# pylint: disable=C0103,R0903,W1514,C0301,E1121,R0913,R0917

""" module docstring """

import logging

from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import lru_cache

from intervaltree import IntervalTree

from . import get_database
from .models import db
from .importers import GqDatabaseImporter


logger = logging.getLogger(__name__)


@dataclass
class OverlapResultHeader:
    has_target: bool


@dataclass
class OverlapTarget:
    start: int = None
    end: int = None
    seq_feature: db.AnnotatedSequence = None
    cov_start: int = None
    cov_end: int = None

    def has_annotation(self):
        return self.seq_feature is not None and bool(self.seq_feature.annotation_str)


class AnnotationDatabaseManager(ABC):
    def __init__(self):
        ...

    @classmethod
    def from_db(cls, db_path, in_memory=True):
        if isinstance(db_path, str):
            return SQL_ADM(get_database(db_path, in_memory=in_memory)[1])
        if isinstance(db_path, GqDatabaseImporter):
            return Dict_ADM(db_path)
        return SQL_ADM(db_path)

    @abstractmethod
    def query_sequence_internal(self, seqid, start=None, end=None):
        ...

    @abstractmethod
    def get_features(self, category_id):
        ...

    @abstractmethod
    def get_categories(self):
        ...

    def query_sequence(self, seqid, start=None, end=None):
        """ Returns strand, seq-feature id (e.g. contig id), functional categories """
        db_sequence = self.query_sequence_internal(seqid, start=start, end=end)

        if db_sequence is None:
            return None
        categories = tuple()
        for cat_features in db_sequence.annotation_str.strip().split(";"):
            category, features = cat_features.split("=")
            categories += ((category, tuple(feature.strip() for feature in features.split(",") if feature.strip())),)
        return db_sequence.strand, db_sequence.featureid, categories

    @abstractmethod
    def query_feature(self, feature_id):
        ...

    @abstractmethod
    def query_category(self, category_id):
        ...

    @lru_cache(maxsize=10000)
    def get_interval_tree(self, seqid):
        db_sequences = self.get_db_sequence(seqid)
        return IntervalTree.from_tuples(
            sorted((seq.start - 1, seq.end) for seq in db_sequences)
        )

    @lru_cache(maxsize=10000)
    @abstractmethod
    def get_db_sequence(self, seqid):
        """ abstract method for sequence retrieval """
        # pylint: disable=W2301
        ...

    def get_interval_tree_overlaps(self, seqid, qstart, qend, with_coverage=False):
        db_sequences = {
            (seq.start, seq.end + 1): seq
            for seq in self.get_db_sequence(seqid)
        }

        # yield bool(db_sequences), None, None, None
        yield OverlapResultHeader(bool(db_sequences))

        itree = IntervalTree.from_tuples(db_sequences)

        for interval in itree[qstart:qend + 1]:
            # yield None, interval.begin, interval.end - 1, db_sequences.get((interval.begin, interval.end))
            cov_start, cov_end = (max(qstart, interval.begin), min(qend, interval.end - 1)) if with_coverage else (None, None)
            yield OverlapTarget(interval.begin, interval.end - 1, db_sequences.get((interval.begin, interval.end)), cov_start, cov_end)

    def get_interval_overlaps(self, seqid, qstart, qend):
        """ return all intervals overlapping the query read """
        db_sequences = self.get_db_sequence(seqid)

        # yield bool(db_sequences), None, None, None
        yield OverlapResultHeader(bool(db_sequences))

        for seq in db_sequences:
            # we're assuming
            # 1) database coordinates in 1-based, closed interval
            # 2) read coordinates in 0-based, closed interval (pysam!)
            if qend < seq.start - 1 or seq.end - 1 < qstart:
                continue

            # yield None, seq.start, seq.end, None
            yield OverlapTarget(seq.start, seq.end, None)

            # interval = seq.start, seq.end
            # sstart, send = seq.start - 1, seq.end

            # if qend < sstart or send < qstart:
            #     continue

            # yield interval

            # if sstart <= qstart <= qend <= send:
            #     yield interval, (qstart, qend)
            # elif qstart < sstart:
            #     yield interval, (sstart, min(qend, send))
            # elif send < qend:
            #     yield interval, (max(qstart, sstart), send)

    @staticmethod
    def calc_covered_fraction(qstart, qend, sstart, send):
        if sstart <= qstart <= qend <= send:
            return (qstart, qend)
        if qstart < sstart:
            return (sstart, min(qend, send))
        if send < qend:
            return (max(qstart, sstart), send)
        raise ValueError(f"Cannot happen. interval=({sstart}, {send}) vs ({qstart}, {qend})")

    def get_overlaps(self, seqid, start, end, domain_mode=False, with_coverage=False):

        # def calc_covered_fraction(start, end, interval):
        #     if interval.begin <= start <= end <= interval.end:
        #         return start, end
        #     if start < interval.begin:
        #         return interval.begin, min(end, interval.end)
        #     if interval.end < end:
        #         return max(start, interval.begin), interval.end
        #     raise ValueError(f"Cannot happen. interval=({interval.begin}, {interval.end}) vs ({start}, {end})")

        if domain_mode:
            return self.get_interval_overlaps(seqid, start, end)

        return self.get_interval_tree_overlaps(seqid, start, end, with_coverage=with_coverage)

        # return (
        #     (interval.begin + 1, interval.end)
        #     for interval in self.get_interval_tree(seqid)[start:end]
        # )

        #     overlaps = self.get_interval_tree(seqid)[start:end]
        #     covered = (
        #         calc_covered_fraction(start, end, interval)
        #         for interval in overlaps
        #     )
        # return overlaps, covered

    # pylint: disable=E1120
    def clear_caches(self):
        logger.info("%s", self.get_interval_tree.cache_info())
        self.get_interval_tree.cache_clear()
        logger.info("%s", self.get_db_sequence.cache_info())
        self.get_db_sequence.cache_clear()


class SQL_ADM(AnnotationDatabaseManager):
    def __init__(self, db_path):
        super().__init__()
        self.db_session = db_path

    def query_sequence_internal(self, seqid, start=None, end=None):
        if start is not None and end is not None:
            db_sequence = self.db_session.query(db.AnnotatedSequence) \
                .filter(db.AnnotatedSequence.seqid == seqid) \
                .filter((db.AnnotatedSequence.start == start) & (db.AnnotatedSequence.end == end)).one_or_none()
        else:
            db_sequence = self.db_session.query(db.AnnotatedSequence).filter(db.AnnotatedSequence.seqid == seqid).one_or_none()

        return db_sequence

    def query_feature(self, feature_id):
        return self.db_session.query(db.Feature).filter(db.Feature.id == feature_id).join(db.Category, db.Feature.category_id == db.Category.id).one_or_none()

    def get_features(self, category_id):
        return self.db_session.query(db.Feature).filter(db.Feature.category_id == category_id).order_by(db.Feature.name).all()

    def query_category(self, category_id):
        return self.db_session.query(db.Category).filter(db.Category.id == category_id).one_or_none()

    def get_categories(self):
        return self.db_session.query(db.Category).all()

    @lru_cache(maxsize=10000)
    def get_db_sequence(self, seqid, start=None, end=None):
        seqs = self.db_session.query(db.AnnotatedSequence).filter(db.AnnotatedSequence.seqid == seqid).all()
        return [
            seq
            for seq in seqs
            if (start is None or start == seq.start) and (end is None or end == seq.end)
        ]


class Dict_ADM(AnnotationDatabaseManager):
    def __init__(self, db_path):
        super().__init__()
        self.db = db_path

    def query_sequence_internal(self, seqid, start=None, end=None):
        if start is not None and end is not None:
            seqs = [
                seq
                for seq in self.db.annotations.get(seqid, [])
                if seq.start == start and seq.end == end
            ] + [None]
            return seqs[0]
        return self.db.annotations.get(seqid, [None])[0]

    def query_feature(self, feature_id):
        return self.db.features.get(int(feature_id))

    def get_features(self, category_id):
        yield from sorted(
            (
                feature
                for feature in self.db.features.values()
                if feature.category_id == category_id
            ),
            key=lambda f: f.name,
        )

    def query_category(self, category_id):
        return self.db.categories.get(int(category_id))

    def get_categories(self):
        yield from self.db.categories.values()

    @lru_cache(maxsize=10000)
    def get_db_sequence(self, seqid, start=None, end=None):
        seqs = self.db.annotations.get(seqid, [])
        return [
            seq
            for seq in seqs
            if (start is None or start == seq.start) and (end is None or end == seq.end)
        ]

    def dump(self, prefix):
        # pylint: disable=C0415
        import pickle
        with open(f"{prefix}.annotations.dat", "wb") as _out:
            pickle.dump(self.db.annotations, _out)
        with open(f"{prefix}.features.dat", "wb") as _out:
            pickle.dump(self.db.features, _out)
        with open(f"{prefix}.categories.dat", "wb") as _out:
            pickle.dump(self.db.categories, _out)

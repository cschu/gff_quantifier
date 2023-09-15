# pylint: disable=C0103,R0903,W1514,C0301,E1121,R0913

""" module docstring """

import logging

from abc import ABC, abstractmethod
from functools import lru_cache

from intervaltree import IntervalTree

from . import get_database
from .models import db
from .db_import import GqDatabaseImporter


logger = logging.getLogger(__name__)


class AnnotationDatabaseManager(ABC):
    def __init__(self):
        ...

    @classmethod
    def from_db(cls, db_path):
        if isinstance(db_path, str):
            return SQL_ADM(get_database(db_path)[1])
        if isinstance(db_path, GqDatabaseImporter):
            return Dict_ADM(db_path)
        return SQL_ADM(db_path)

    @abstractmethod
    def query_sequence_internal(self, seqid, start=None, end=None):
        ...

    def query_sequence(self, seqid, start=None, end=None):
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

    @abstractmethod
    def get_db_sequence(self, seqid):
        """ abstract method for sequence retrieval """
        # pylint: disable=W2301
        ...

    def get_interval_overlaps(self, seqid, qstart, qend):
        """ return all intervals overlapping the query read """
        db_sequences = self.get_db_sequence(seqid)

        for seq in db_sequences:
            # we're assuming
            # 1) database coordinates in 1-based, closed interval
            # 2) read coordinates in 0-based, closed interval (pysam!)
            if qend < seq.start - 1 or seq.end - 1 < qstart:
                continue

            yield seq.start, seq.end

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

    def get_overlaps(self, seqid, start, end, domain_mode=False):

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
        return (
            (interval.begin + 1, interval.end)
            for interval in self.get_interval_tree(seqid)[start:end]
        )
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

    def query_category(self, category_id):
        return self.db_session.query(db.Category).filter(db.Category.id == category_id).one_or_none()

    @lru_cache(maxsize=10000)
    def get_db_sequence(self, seqid):
        return self.db_session.query(db.AnnotatedSequence).filter(db.AnnotatedSequence.seqid == seqid).all()


class Dict_ADM(AnnotationDatabaseManager):
    def __init__(self, db_path):
        super().__init__()
        self.db = db_path

    def query_sequence_internal(self, seqid, start=None, end=None):
        if start is not None and end is not None:
            seqs = [seq for seq in self.db.annotations.get(seqid, [None]) if seq.start == start and seq.end == end]
            return seqs[0]
        return self.db.annotations.get(seqid, [None])[0]

    def query_feature(self, feature_id):
        return self.db.features.get(int(feature_id))

    def query_category(self, category_id):
        return self.db.categories.get(int(category_id))

    @lru_cache(maxsize=10000)
    def get_db_sequence(self, seqid):
        return self.db.annotations.get(seqid, [])

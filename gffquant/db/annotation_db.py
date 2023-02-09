# pylint: disable=C0103,R0903,W1514,C0301,E1121,R0913

""" module docstring """

import logging

from functools import lru_cache

from intervaltree import IntervalTree

from gffquant.db import get_database
from gffquant.db.models import db

logger = logging.getLogger(__name__)


class AnnotationDatabaseManager:
    def __init__(self, db_path):
        _, self.dbsession = get_database(db_path)

    def query_sequence(self, seqid, start=None, end=None):
        if start is not None and end is not None:
            db_sequence = self.dbsession.query(db.AnnotatedSequence) \
                .filter(db.AnnotatedSequence.seqid == seqid) \
                .filter((db.AnnotatedSequence.start == start) & (db.AnnotatedSequence.end == end)).one_or_none()
        else:
            db_sequence = self.dbsession.query(db.AnnotatedSequence).filter(db.AnnotatedSequence.seqid == seqid).one_or_none()
        if db_sequence is None:
            return None
        categories = tuple()
        for cat_features in db_sequence.annotation_str.strip().split(";"):
            category, features = cat_features.split("=")
            categories += ((category, tuple(feature.strip() for feature in features.split(",") if feature.strip())),)
        return db_sequence.strand, db_sequence.featureid, categories

    def query_feature(self, feature_id):
        return self.dbsession.query(db.Feature).filter(db.Feature.id == feature_id).join(db.Category, db.Feature.category_id == db.Category.id).one_or_none()

    def query_category(self, category_id):
        return self.dbsession.query(db.Category).filter(db.Category.id == category_id).one_or_none()

    @lru_cache(maxsize=10000)
    def get_interval_tree(self, seqid):
        db_sequences = self.get_db_sequence(seqid)
        return IntervalTree.from_tuples(
            sorted((seq.start - 1, seq.end) for seq in db_sequences)
        )

    @lru_cache(maxsize=10000)
    def get_db_sequence(self, seqid):
        return self.dbsession.query(db.AnnotatedSequence).filter(db.AnnotatedSequence.seqid == seqid).all()

    def get_interval_overlaps(self, seqid, qstart, qend, calc_coverage=True):
        db_sequences = self.get_db_sequence(seqid)

        for seq in db_sequences:
            interval = seq.start, seq.end
            sstart, send = seq.start - 1, seq.end

            if qend < sstart or send < qstart:
                continue
            
            covered_interval = self.calc_covered_fraction(qstart, qend, sstart, send) if calc_coverage else interval
            yield interval, covered_interval
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

    def get_overlaps(self, seqid, start, end, domain_mode=False, calc_coverage=False):

        # def calc_covered_fraction(start, end, interval):
        #     if interval.begin <= start <= end <= interval.end:
        #         return start, end
        #     if start < interval.begin:
        #         return interval.begin, min(end, interval.end)
        #     if interval.end < end:
        #         return max(start, interval.begin), interval.end
        #     raise ValueError(f"Cannot happen. interval=({interval.begin}, {interval.end}) vs ({start}, {end})")

        if domain_mode:
            return self.get_interval_overlaps(seqid, start, end, calc_coverage=calc_coverage)
        return (
            (
                (interval.begin + 1, interval.end),
                (self.calc_covered_fraction(start, end, interval.begin, interval.end) if calc_coverage else (interval.begin + 1, interval.end))
            )
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

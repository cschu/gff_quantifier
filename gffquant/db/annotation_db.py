# pylint: disable=C0103,R0903,W1514,C0301

""" module docstring """

from functools import lru_cache

from intervaltree import IntervalTree

from gffquant.db import get_database
from gffquant.db.models import db


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
        db_sequences = self.dbsession.query(db.AnnotatedSequence).filter(db.AnnotatedSequence.seqid == seqid).all()
        return IntervalTree.from_tuples(
            sorted((seq.start, seq.end) for seq in db_sequences)
        )

    def get_overlaps(self, seqid, start, end):

        def calc_covered_fraction(start, end, interval):
            if interval.begin <= start <= end <= interval.end:
                return start, end
            if start < interval.begin:
                return interval.begin, min(end, interval.end)
            if interval.end < end:
                return max(start, interval.begin), interval.end
            raise ValueError(f"Cannot happen. interval=({interval.begin}, {interval.end}) vs ({start}, {end})")

        overlaps = self.get_interval_tree(seqid)[start:end]
        covered = (
            calc_covered_fraction(start, end, interval)
            for interval in overlaps
        )
        return overlaps, covered


"""
    def get_overlaps(self, ref, start, end, cache_data=False):
        def calc_covered_fraction(start, end, interval):
            cstart, cend = start, end
            if interval.begin <= start <= end <= interval.end:
                ...
            elif start < interval.begin:
                cstart, cend = interval.begin, min(end, interval.end)
            elif interval.end < end:
                cstart, cend = max(start, interval.begin), interval.end
            return cstart, cend

        overlaps = self._get_tree(ref, cache_data=cache_data)[start:end]
        covered = [calc_covered_fraction(start, end, interval) for interval in overlaps]

        return overlaps, covered

    @lru_cache(maxsize=4096)
    def _get_tree(self, ref, cache_data=False):
        if self.reference_type == "domain":
            return IntervalTree.from_tuples(
                sorted((start, end) for start, end in self.db.get(ref, {}))
            )

        return IntervalTree.from_tuples(
            sorted(
                [key[1:] for key in self._read_data(ref, include_payload=cache_data)]
            )
        )




@lru_cache(maxsize=4096)
    def _read_data(self, ref_id, include_payload=False):
        gff_annotation = {}
        for offset, size in self.db_index.get(ref_id, []):
            self.db.seek(offset)
            for line in self.db.read(size).strip("\n").split("\n"):
                if not line.startswith("#"):
                    line = line.strip().split("\t")
                    features = {}
                    if include_payload:
                        features = (("strand", line[6]),)
                        features += tuple(
                            (
                                item.split("=")[0],
                                tuple(sorted(item.split("=")[1].split(","))),
                            )
                            for item in line[8].strip().split(";")
                        )
                    key = (line[0], int(line[3]), int(line[4]) + 1)
                    gff_annotation[key] = features
        if not gff_annotation and not include_payload:
            print(
                f"WARNING: contig {ref_id} does not have an annotation in the index.",
                file=sys.stderr,
                flush=True,
            )
        return gff_annotation

    def get_data(self, ref, start, end):
        if self.reference_type == "domain":
            dom_features = self.db.get(ref, {}).get((start, end), [])
            features = (("strand", None), ("ID", ref))
            features += (("domtype", tuple(dom_features)),)
            return features if dom_features else {}

        return self._read_data(ref, include_payload=True).get((ref, start, end), {})

    # pylint: disable=E1120
    def clear_caches(self):
        print(self._read_data.cache_info(), flush=True)
        self._read_data.cache_clear()
        print(self._get_tree.cache_info(), flush=True)
        self._get_tree.cache_clear()

"""

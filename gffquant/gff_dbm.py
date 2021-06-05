import sys
import gzip
from functools import lru_cache

from intervaltree import IntervalTree


class GffDatabaseManager:

    def iterate(self):
        header = None
        with self.db as db_stream:
            for line in db_stream:
                line = line.strip()
                if line.startswith("#"):
                    header = line.strip("#").split("\t")
                else:
                    line = line.split("\t")
                    features = list()
                    for feat_cat, subfeatures in zip(header[6:], line[6:]):
                        subfeatures = tuple(sf for sf in subfeatures.strip().split(",") if sf)
                        if subfeatures:
                            features.append((feat_cat, subfeatures))
                    yield line[0], (("strand", None),) + tuple(features)

    def _read_index(self, f):
        db_index = dict()
        for line in open(f, "rt"):
            line = line.strip().split("\t")
            db_index.setdefault(line[0], list()).append(list(map(int, line[1:3])))
        return db_index

    def __init__(self, db, db_index=None):
        gz_magic = b"\x1f\x8b\x08"
        gzipped = open(db, "rb").read(3).startswith(gz_magic)
        if db_index:
            if gzipped:
                raise ValueError(f"Database {db} is gzipped. This doesn't work together with an index. Please unzip and re-index.")
            _open = open
            self.db_index = self._read_index(db_index)
        else:
            _open = gzip.open if gzipped else open
            self.db_index = None
        self.db = _open(db, "rt")

    @lru_cache(maxsize=4096)
    def _read_data(self, ref_id, include_payload=False):
        gff_annotation = dict()
        for offset, size in self.db_index.get(ref_id, list()):
            self.db.seek(offset)
            for line in self.db.read(size).strip("\n").split("\n"):
                if not line.startswith("#"):
                    line = line.strip().split("\t")
                    features = dict()
                    if include_payload:
                        features = (("strand", line[6]),)
                        features += tuple((item.split("=")[0], tuple(sorted(item.split("=")[1].split(",")))) for item in line[8].strip().split(";"))
                    key = (line[0], int(line[3]), int(line[4]) + 1)
                    gff_annotation[key] = features
        if not gff_annotation and not include_payload:
            print("WARNING: contig {contig} does not have an annotation in the index.".format(contig=ref_id), file=sys.stderr, flush=True)
        return gff_annotation

    @lru_cache(maxsize=4096)
    def _get_tree(self, ref, cache_data=False):
        return IntervalTree.from_tuples(sorted([key[1:] for key in self._read_data(ref, include_payload=cache_data)]))

    def get_data(self, ref, start, end):
        return self._read_data(ref, include_payload=True).get((ref, start, end), dict())

    def get_overlaps(self, ref, start, end, cache_data=False):
        return self._get_tree(ref, cache_data=cache_data)[start:end]

    def clear_caches(self):
        print(self._read_data.cache_info(), flush=True)
        self._read_data.cache_clear()
        print(self._get_tree.cache_info(), flush=True)
        self._get_tree.cache_clear()

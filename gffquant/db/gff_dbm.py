# pylint: disable=C0103,R0903,W1514,R1702

""" module docstring """

import sys
import gzip
from functools import lru_cache

from intervaltree import IntervalTree


class EmapperFormat:
    def __init__(self, query_field, fields, categories, na):
        self.query_field = query_field
        self.categories = dict(zip(fields, categories))
        self.na = na

    def get_category(self, index):
        return self.categories.get(index)

    def parse_emapper_line(self, line):
        line = line.strip().split("\t")
        categories = []
        for i, col in enumerate(line):
            category = self.get_category(i)
            if category is not None:
                col = col.strip()
                if col and col != self.na:
                    features = tuple(
                        feature.strip()
                        for feature in col.split(",")
                        if feature.strip()
                    )
                    if features:
                        categories.append((category, features))
        if categories:
            return line[self.query_field], (("strand", None),) + tuple(categories)

        return None


EMAPPER_FORMATS = {
    "v1": EmapperFormat(
        0,
        (5, 6, 7, 9, 11),
        (
            "Gene_Ontology_terms",
            "KEGG_ko",
            "BiGG_Reaction",
            "eggNOG_OGs",
            "COG_Functional_Category",
        ),
        "",
    ),
    "v2": EmapperFormat(
        0,
        (6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20),
        (
            "Gene_Ontology_terms",
            "EC_number",
            "KEGG_ko",
            "KEGG_Pathway",
            "KEGG_Module",
            "KEGG_Reaction",
            "KEGG_rclass",
            "BRITE",
            "KEGG_TC",
            "CAZy",
            "BiGG_Reaction",
            "eggNOG_OGs",
            "COG_Functional_Category",
        ),
        "",
    ),
    "v2.1.2": EmapperFormat(
        0,
        (4, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,),
        (
            "eggNOG_OGs",
            "COG_cat",
            "GOs",
            "EC",
            "KEGG_ko",
            "KEGG_Pathway",
            "KEGG_Module",
            "KEGG_Reaction",
            "KEGG_rclass",
            "BRITE",
            "KEGG_TC",
            "CAZy",
            "BiGG_Reaction",
            "PFAMs",
        ),
        "-",
    ),
}


class GffDatabaseManager:
    def iterate(self, bufsize=4000000):
        with self.db as db_stream:
            tail = ""
            while True:
                chunk = "".join((tail, db_stream.read(bufsize).decode()))
                if not chunk:
                    break
                chunk = chunk.split("\n")
                chunk, tail = chunk[:-1], chunk[-1]
                for line in chunk:
                    if line[0] != "#":
                        annotation = self.emapper_format.parse_emapper_line(line)
                        if annotation is not None:
                            yield annotation
                        # line = line.strip().split("\t")
                        # categories = []
                        # for i, col in enumerate(line):
                        #     category = self.emapper_format.get_category(i)
                        #     if category is not None:
                        #         col = col.strip().split(",")
                        #         if col != self.emapper_format.na:
                        #             features = tuple(
                        #                 feature.strip()
                        #                 for feature in col
                        #                 if feature.strip()
                        #             )
                        #             if features:
                        #                 categories.append((category, features))
                        # if categories:
                        #     yield line[self.emapper_format.query_field], (
                        #         ("strand", None),
                        #     ) + tuple(categories)

    def _read_index(self, f):
        self.db_index = {}
        with open(f, "rt") as _in:
            for line in _in:
                line = line.strip().split("\t")
                self.db_index.setdefault(line[0], []).append(list(map(int, line[1:3])))

    def __init__(self, db, reference_type, db_index=None, emapper_version="v2"):
        gz_magic = b"\x1f\x8b\x08"
        # pylint: disable=R1732,W0511
        gzipped = open(db, "rb").read(3).startswith(gz_magic)
        # TODO: can dbm be written as contextmanager?
        _open = gzip.open if gzipped else open
        self.reference_type = reference_type
        if db_index:
            if gzipped:
                raise ValueError(
                    f"Database {db} is gzipped. "
                    "This doesn't work together with an index. Please unzip and re-index."
                )
            self._read_index(db_index)
            self.db = _open(db, "rt")
        elif self.reference_type == "domain":  # bed
            self.db = {}
            for line in _open(db, "rt"):
                line = line.strip().split("\t")
                self.db.setdefault(line[0], {}).setdefault(
                    (int(line[1]), int(line[2])), []
                ).append(line[3])
        else:
            _open = gzip.open if gzipped else open
            self.db = _open(db, "rb")
            self.db_index = None

        self.emapper_format = EMAPPER_FORMATS.get(emapper_version)
        if self.reference_type in ("gene", "genes") and not self.emapper_format:
            raise ValueError(
                f"Cannot find emapper parse instructions for version {emapper_version}."
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

    def get_data(self, ref, start, end):
        if self.reference_type == "domain":
            dom_features = self.db.get(ref, {}).get((start, end), [])
            features = (("strand", None), ("ID", ref))
            features += (("domtype", tuple(dom_features)),)
            return features if dom_features else {}

        return self._read_data(ref, include_payload=True).get((ref, start, end), {})

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

    # pylint: disable=E1120
    def clear_caches(self):
        print(self._read_data.cache_info(), flush=True)
        self._read_data.cache_clear()
        print(self._get_tree.cache_info(), flush=True)
        self._get_tree.cache_clear()

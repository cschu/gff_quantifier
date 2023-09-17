# pylint: disable=C0103,R0902,R0913

""" module docstring """

import contextlib
import gzip
import json

from abc import ABC, abstractmethod

from .models import db


class GqDatabaseImporter(ABC):
    def __init__(self, logger, input_data, db_path=None, db_session=None):
        self.logger = logger
        self.db_path = db_path
        self.db_session = db_session
        self.code_map = {}
        self.annotations = {}
        self.nseqs = 0
        self.categories = {}
        self.features = {}

        self.gather_category_and_feature_data(input_data)
        self.process_annotations(input_data)

    @staticmethod
    def get_open_function(f):
        gz_magic = b"\x1f\x8b\x08"
        # pylint: disable=R1732,W0511
        gzipped = open(f, "rb").read(3).startswith(gz_magic)
        return gzip.open if gzipped else open

    @abstractmethod
    def parse_categories(self, input_data):
        ...

    @abstractmethod
    def parse_annotations(self, input_data):
        ...

    def gather_category_and_feature_data(self, input_data):
        self.logger.info("First pass: gathering category and feature information.")

        _open = GqDatabaseImporter.get_open_function(input_data)

        with _open(input_data, "rt") as _in:
            cat_d = self.parse_categories(_in)

        self.logger.info("Building code map and dumping category and feature encodings.")

        feature_offset = 0
        _map_out = gzip.open(self.db_path + ".code_map.json.gz", "wt") if self.db_path else contextlib.nullcontext()
        with _map_out:
            for category, features in sorted(cat_d.items()):
                self.code_map[category] = {
                    "key": len(self.code_map),
                    "features": {
                        feature: (i + feature_offset) for i, feature in enumerate(sorted(features))
                    }
                }
                feature_offset += len(features)

                db_category = db.Category(id=self.code_map[category]["key"], name=category)
                self.categories[db_category.id] = db_category
                if self.db_session is not None:
                    self.db_session.add(db_category)
                    self.db_session.commit()

                for feature, fid in self.code_map[category]["features"].items():
                    db_feature = db.Feature(id=fid, name=feature, category=db_category)
                    self.features[db_feature.id] = db_feature
                    if self.db_session is not None:
                        self.db_session.add(db_feature)

                if self.db_session is not None:
                    self.db_session.commit()

            if self.db_path:
                json.dump(self.code_map, _map_out)

    def process_annotations(self, input_data):
        self.logger.info("Second pass: Encoding sequence annotations")

        _open = GqDatabaseImporter.get_open_function(input_data)

        single_category = getattr(self, "single_category") if hasattr(self, "single_category") else "domain"

        with _open(input_data, "rt") as _in:
            annotation_data = self.parse_annotations(_in)

            i = 0
            for i, ((gid, start, end), features) in enumerate(annotation_data.items(), start=1):
                if i % 100000 == 0:
                    if self.nseqs:
                        self.logger.info("    Loaded %s entries. (%s%%)", i, round(i / self.nseqs * 100, 3))
                    else:
                        self.logger.info("    Loaded %s entries.", str(i))

                encoded = []
                enc_category = self.code_map[single_category]['key']
                enc_features = sorted(self.code_map[single_category]['features'][feature] for feature in features)
                encoded.append((enc_category, ",".join(map(str, enc_features))))
                encoded = ";".join(f"{cat}={features}" for cat, features in sorted(encoded))

                db_sequence = db.AnnotatedSequence(
                    seqid=gid,
                    featureid=None,
                    start=int(start),
                    end=int(end),
                    annotation_str=encoded,
                )

                self.annotations.setdefault(gid, []).append(db_sequence)

                if self.db_session:
                    self.db_session.add(db_sequence)

            if self.nseqs:
                self.logger.info("    Loaded %s entries. (%s%%)", i, round(i / self.nseqs * 100, 3))
            else:
                self.logger.info("    Loaded %s entries.", str(i))

            if self.db_session:
                self.db_session.commit()


# class DomainBedDatabaseImporter(GqDatabaseImporter):
class SmallDatabaseImporter(GqDatabaseImporter):
    def __init__(
        self,
        logger,
        input_data,
        db_path=None,
        db_session=None,
        single_category="domain",
        sep="\t",
        coords="bed",
    ):
        self.single_category = single_category
        self.sep = sep
        # we store everything as 1-based, closed intervals internally
        # bed coords coming in as [x,y)_0 -> [x+1, y]_1
        self.coordinate_modifiers = (1, 0) if coords == "bed" else (0, 0)
        self.cols = (0,1,2,3) if coords == "bed" else (0,1,2,4)
        # GMGC10.000_000_128.UNKNOWN,1,420,4.958096936477401e-50,GH88,417

        super().__init__(logger, input_data, db_path=db_path, db_session=db_session)

    def parse_categories(self, _in):
        categories = {}

        for self.nseqs, line in enumerate(_in, start=1):
            line = line.strip().split(self.sep)
            categories.setdefault(self.single_category, set()).update(line[self.cols[-1]].split(","))

        self.logger.info("    Parsed %s entries.", self.nseqs)
        return categories

    def parse_annotations(self, _in):
        annotations = {}
        for i, line in enumerate(_in, start=1):
            if i % 10000 == 0 and self.db_session:
                self.db_session.commit()
            line = line.strip().split(self.sep)
            gid, start, end, features = (c for i, c in enumerate(line) if i in self.cols)
            # we store everything as 1-based, closed intervals internally
            start = int(start) + self.coordinate_modifiers[0]
            end = int(end) + self.coordinate_modifiers[1]
            # was: start + 1
            annotations.setdefault((gid, start, end), set()).update(features.split(","))

        return annotations

# pylint: disable=C0103,R0902,R0913,W2301

""" module docstring """

import gzip
import logging

from abc import ABC, abstractmethod

from ..models import db


logger = logging.getLogger(__name__)


class GqDatabaseImporter(ABC):
    """ Base database importer class"""
    def __init__(self, input_data, db_path=None, db_session=None, na_char="-"):
        self.db_path = db_path
        self.db_session = db_session
        self.code_map = {}
        self.annotations = {}
        self.nseqs = 0
        self.categories = {}
        self.features = {}
        self.open_function = GqDatabaseImporter.get_open_function(input_data)
        self.na_char = na_char

        # self.gather_category_and_feature_data(input_data)
        # self.process_annotations(input_data)
        self.build_database(input_data)

    @staticmethod
    def get_open_function(f):
        """ Returns a file open function corresponding to gzip-compression status. """
        gz_magic = b"\x1f\x8b\x08"
        # pylint: disable=R1732,W0511
        gzipped = open(f, "rb").read(3).startswith(gz_magic)
        return gzip.open if gzipped else open

    @abstractmethod
    def parse_annotations(self, input_data):
        """ abstract method to parse annotations from various data formats """
        ...

    def build_database(self, input_data):

        category_map, feature_map = {}, {}

        with self.open_function(input_data, "rb") as _in:
            annotation_data = self.parse_annotations(_in)

            i = 0
            for i, (seq_feature, annotation) in enumerate(annotation_data, start=1):
                if i % 100000 == 0:
                    if self.db_session is not None:
                        self.db_session.commit()

                    logger.info("    Loaded %s entries.", str(i))

                encoded = []
                for category, features in annotation:
                    cat_enc = category_map.setdefault(category, len(category_map))

                    feat_enc = sorted(
                        feature_map.setdefault((cat_enc, feat), len(feature_map))
                        for feat in features
                    )

                    encoded.append((cat_enc, ",".join(map(str, feat_enc))))

                seq_feature.annotation_str = ";".join(
                    f"{cat}={features}" for cat, features in sorted(encoded)
                )

                if self.db_session:
                    self.db_session.add(seq_feature)
                else:
                    self.annotations.setdefault(seq_feature.seqid, []).append(seq_feature)

            logger.info("Finished loading %s entries.", str(i))

            logger.info("Loading categories and features.")

            self.categories.update(
                {
                    cat_id: db.Category(id=cat_id, name=cat_name)
                    for cat_name, cat_id in category_map.items()
                }
            )

            self.features.update(
                {
                    feat_id: db.Feature(id=feat_id, name=feat_name, category_id=cat_id)
                    for (cat_id, feat_name), feat_id in feature_map.items()
                }
            )

            if self.db_session is not None:
                for category in self.categories.values():
                    self.db_session.add(category)
                for feature in self.features.values():
                    self.db_session.add(feature)

                self.db_session.commit()

            logger.info("Finished loading categories and features.")

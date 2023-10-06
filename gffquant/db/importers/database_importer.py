# pylint: disable=C0103,R0902,R0913,W2301,R0914

""" module docstring """

import gzip
import logging

from abc import ABC, abstractmethod

from sqlalchemy import insert

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
        annotations = []

        with self.open_function(input_data, "rb") as _in:
            annotation_data = self.parse_annotations(_in)

            i = 0
            for i, (seq_feature, annotation) in enumerate(annotation_data, start=1):
                if not annotation:
                    continue

                if i % 100000 == 0:
                    if self.db_session is not None and annotations:
                        self.db_session.execute(
                            insert(db.AnnotatedSequence),
                            [ann.__dict__ for ann in annotations],
                        )
                        # self.db_session.bulk_save_objects(annotations)
                        self.db_session.commit()
                        annotations.clear()
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

                annotations.append(seq_feature)

            logger.info("    Loaded %s entries.", str(self.nseqs))

            self.categories = {
                cat_id: db.Category(id=cat_id, name=cat_name)
                for cat_name, cat_id in category_map.items()
            }

            self.features = {
                feat_id: db.Feature(id=feat_id, name=feat_name, category_id=cat_id)
                for (cat_id, feat_name), feat_id in feature_map.items()
            }

            if self.db_session is not None:
                self.db_session.execute(
                    insert(db.AnnotatedSequence),
                    [ann.__dict__ for ann in annotations],
                )
                # bulk_save_objects(annotations)

                self.db_session.execute(
                    insert(db.Category),
                    [cat.__dict__ for cat in self.categories.values()],
                )
                # self.db_session.bulk_save_objects(self.categories.values())

                self.db_session.execute(
                    insert(db.Feature),
                    [feat.__dict__ for feat in self.features.values()],
                )
                # self.db_session.bulk_save_objects(self.features.values())
                self.db_session.commit()

                annotations.clear()
                self.features.clear()
                self.categories.clear()
            else:
                for seq_feature in annotations:
                    self.annotations.setdefault(seq_feature.seqid, []).append(seq_feature)

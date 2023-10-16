# pylint: disable=C0103,R0902,R0913,W2301

""" module docstring """

import logging

from .database_importer import GqDatabaseImporter
from ..models import db
from .settings import DB_SETTINGS_SELECTION


logger = logging.getLogger(__name__)


class SmallDatabaseImporter(GqDatabaseImporter):
    """ Importer for small dict-based databases. """
    def __init__(
        self,
        db_path=None,
        db_session=None,
        single_category="domain",
        db_format="default",
    ):
        self.single_category = single_category

        self.db_settings = DB_SETTINGS_SELECTION.get(db_format)
        if self.db_settings is None:
            raise ValueError(f"{db_format=} is not recognised.")

        super().__init__(db_path=db_path, db_session=db_session)

    def parse_annotations(self, input_data, input_data2=None):
        for line in input_data:
            line = line.decode().strip().split(self.db_settings.separator)
            gid, start, end, features = (c for i, c in enumerate(line) if i in self.db_settings.columns)

            # we store everything as 1-based, closed intervals internally
            start = int(start) + self.db_settings.offsets[0]
            end = int(end) + self.db_settings.offsets[1]

            seq_feature = db.AnnotatedSequence(
                seqid=gid,
                featureid=None,
                start=start,
                end=end,
            )

            annotation = [
                (self.single_category, set(features.split(",")))
            ]

            yield seq_feature, annotation

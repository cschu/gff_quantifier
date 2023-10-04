# pylint: disable=C0103,R0902,R0913,W2301

""" module docstring """

import logging

from .database_importer import GqDatabaseImporter
from ..gff_dbm import GffDatabaseManager
from ..models import db


logger = logging.getLogger(__name__)


class GqGeneDatabaseImporter(GqDatabaseImporter):
    def __init__(
        self,
        input_data,
        db_path=None,
        db_session=None,
        emapper_version=None,
        dbm_buffersize=4000000000
    ):
        if emapper_version is None:
            raise ValueError("Missing emapper_version.")
        self.emapper_version = emapper_version
        self.dbm_buffersize = dbm_buffersize

        super().__init__(input_data, db_path=db_path, db_session=db_session)

    def parse_categories(self, input_data):
        categories = {}
        gffdbm = GffDatabaseManager(
            input_data,
            "genes",
            emapper_version=self.emapper_version
        )

        for self.nseqs, (_, region_annotation) in enumerate(
            gffdbm.iterate(bufsize=self.dbm_buffersize), start=1
        ):
            for category, features in region_annotation[1:]:
                categories.setdefault(category, set()).update(features)

        return categories

    def parse_annotations(self, input_data):
        gffdbm = GffDatabaseManager(
            input_data,
            "genes",
            emapper_version=self.emapper_version,
        )
        for ref, region_annotation in gffdbm.iterate(bufsize=self.dbm_buffersize):
            _, strand = region_annotation[0]

            seq_feature = db.AnnotatedSequence(
                seqid=ref,
                featureid=None,
                strand=int(strand == "+") if strand is not None else None,
            )
            annotation = {
                category: set(features).difference({"-"})
                for category, features in region_annotation[1:]
            }

            yield seq_feature, annotation

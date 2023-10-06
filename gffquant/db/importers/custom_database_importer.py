# pylint: disable=C0103,R0902,R0913,W2301

""" module docstring """

import logging

from .database_importer import GqDatabaseImporter
from ..models import db


logger = logging.getLogger(__name__)


class GqCustomDatabaseImporter(GqDatabaseImporter):
    def __init__(
        self,
        input_data,
        db_path=None,
        db_session=None,
        columns=None,
        header=None,
        delimiter="\t",
    ):
        self.columns = columns
        self.header = header
        self.delimiter = delimiter

        super().__init__(input_data, db_path=db_path, db_session=db_session)

    def parse_annotations(self, _in):
        if self.header:
            _ = [next(_in) for _ in range(self.header - 1)]
        header_line = next(_in).strip().strip("#").split(self.delimiter)
        category_cols = set(self.columns) if self.columns is not None else header_line[1:]
        logging.info("    Got header: %s", header_line)
        logging.info("    Got columns: %s", category_cols)
        for col in category_cols:
            if col not in header_line:
                logging.error("    column %s is not present in headers", col)
                raise ValueError(f"column {col} is not present in headers.")
        for self.nseqs, line in enumerate(_in, start=1):
            if line and line[0] != "#":
                line = line.strip().split(self.delimiter)
                line_d = dict(zip(header_line, line))

                seq_feature = db.AnnotatedSequence(
                    seqid=line[0],
                    featureid=None,
                )

                annotation = (
                    (
                        category,
                        set(features.split(",")).difference({"-"})
                    )
                    for category, features in line_d.items()
                    if category in category_cols
                )
            
                annotation = [
                    (category, features)
                    for category, features in annotation
                    if features
                ]
            
                if annotation:
                    yield seq_feature, annotation

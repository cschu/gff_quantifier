# pylint: disable=C0103,R0902,R0913,W2301,W1203,R0917

""" module docstring """

import logging

from .database_importer import GqDatabaseImporter
from ..models import db


logger = logging.getLogger(__name__)


class GqCustomDatabaseImporter(GqDatabaseImporter):
    def __init__(
        self,
        db_path=None,
        db_session=None,
        columns=None,
        skip_header_lines=0,
        header=None,
        delimiter="\t",
    ):

        self.columns = columns
        self.header = header.split(",") if header is not None else None
        self.delimiter = delimiter
        self.skip_header_lines = skip_header_lines

        super().__init__(db_path=db_path, db_session=db_session)

    def _validate_category_columns(self, header_line):
        category_cols = set(self.columns) if self.columns is not None else header_line[1:]
        logging.info("    Got header: %s", header_line)
        logging.info("    Got columns: %s", category_cols)

        if self.columns is not None:
            for col in category_cols:
                if col not in header_line:
                    msg = f"column {col} is not present in headers."
                    logging.error(msg)
                    raise ValueError(msg)

        return category_cols

    def parse_annotations(self, input_data, input_data2=None):
        if self.skip_header_lines > 0:
            try:
                _ = [next(input_data) for _ in range(self.skip_header_lines - 1)]
            except StopIteration as exc:
                msg = f"Reached end of annotation file while skipping header comments ({self.header})."
                logging.error(f"    {msg}")
                raise ValueError(msg) from exc

        if self.header is None:
            try:
                header_line = next(input_data).decode().strip().strip("#").split(self.delimiter)
            except StopIteration as exc:
                msg = "Reached end of annotation file while parsing header line."
                logging.error(f"    {msg}")
                raise ValueError(msg) from exc
        else:
            header_line = self.header

        category_cols = self._validate_category_columns(header_line)

        for self.nseqs, line in enumerate(input_data, start=1):
            line = line.decode()
            if line and line[0] != "#":
                line = line.strip().split(self.delimiter)
                columns = {
                    colname: coldata.strip()
                    for colname, coldata in zip(header_line, line)
                    if colname in category_cols
                }

                # annotation = [
                #     (category, tuple(features.split(",")))
                #     for category, features in columns.items()
                #     if features and features != self.na_char
                # ]
                annotation = tuple(self.extract_features(columns))

                if annotation:
                    seq_feature = db.AnnotatedSequence(
                        seqid=line[0],
                        featureid=None,
                    )
                    yield seq_feature, annotation

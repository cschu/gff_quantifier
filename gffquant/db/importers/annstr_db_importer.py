# pylint: disable=C0103,R0902,R0913,W2301,W1203

""" module docstring """

import gzip
import hashlib
import logging

from .custom_database_importer import GqCustomDatabaseImporter
from ..models import db


logger = logging.getLogger(__name__)


class AnnstrDatabaseImporter(GqCustomDatabaseImporter):
    def __init__(
        self,
        input_data,
        db_path=None,
        db_session=None,
        columns=None,
        seq_column=None,
        skip_header_lines=0,
        header=None,
        delimiter="\t",
    ):
        self.seq_column = seq_column

        super().__init__(
            input_data,
            db_path=db_path,
            db_session=db_session,
            columns=columns,
            header=header,
            delimiter=delimiter,
            skip_header_lines=skip_header_lines,
        )

    def parse_annotations(self, _in):

        if self.header is None:
            header_line = next(_in, None)
            if header_line is None:
                msg = "Reached end of annotation file while parsing header line."
                logging.error(f"    {msg}")
                raise ValueError(msg)

            header_line = header_line.decode()
        else:
            header_line = self.header

        category_cols = self._validate_category_columns(header_line)

        if self.seq_column is not None:
            if self.seq_column not in header_line:
                msg = f"column `{self.seq_column}` is not present in headers.\n{self.header}"
                logging.error(msg)
                raise ValueError(msg)

        annotation_suffices = {}

        def get_ann_hash(s):
            h = hashlib.sha256()
            h.update(s.encode())
            return h.hexdigest()
        
        empty_ann_hash = get_ann_hash("")
        logger.info(f"NO_ANNOTATION => {empty_ann_hash}")

        with gzip.open(f"{self.db_path.replace('sqlite3', 'ffn.gz')}", "wt") as seq_out:
            for self.nseqs, line in enumerate(_in, start=1):
                line = line.decode()
                line = line.strip().split(self.delimiter)
                line_d = {
                    colname: value.strip()
                    for colname, value in zip(header_line + [self.seq_column], line)
                    if colname in category_cols or colname == self.seq_column
                }

                annotation = tuple(
                    (category, tuple(set(sorted(features.split(",")))))
                    for category, features in line_d.items()
                    if features != self.na_char and features and category != self.seq_column
                )
                if not annotation:
                    logger.info(f"Gene `{line[0]}` has no annotation.")

                ann_str = ";".join(
                    f"{category}={','.join(features)}"
                    for category, features in annotation
                )                

                ann_sfx = annotation_suffices.setdefault(
                    ann_str,
                    get_ann_hash(ann_str)
                )

                print(
                    f">{line[0]}.{ann_sfx}", line_d[self.seq_column],
                    sep="\n",
                    file=seq_out
                )

                annotation_string = db.AnnotationString(
                    annotation_hash=ann_sfx,
                )

                yield annotation_string, annotation

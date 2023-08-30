# pylint: disable=R0914,C0103,R0913,W1514
# pylint: disable=duplicate-code

""" module docstring """

import argparse
import gzip
import json
import logging
import sqlite3

from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

from ..db import initialise_db
from ..db.models import db
from ..db.models.meta import Base


logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s'
)


def get_database(db_path):
    engine = create_engine(f"sqlite:///{db_path}")

    db_session = scoped_session(
        sessionmaker(
            autocommit=False,
            autoflush=False,
            enable_baked_queries=True,
            bind=engine
        )
    )

    Base.query = db_session.query_property()

    return engine, db_session


def gather_category_and_feature_data(
    input_data,
    db_path,
    db_session=None,
    columns=None,
    header=None,
    delimiter="\t",
    is_gzipped=False
):
    logging.info("First pass: gathering category and feature information.")

    cat_d = {}
    n = 0
    with (gzip.open if is_gzipped else open)(input_data, "rt") as _in:
        if header:
            _ = [next(_in) for _ in range(header - 1)]
        header_line = next(_in).strip().strip("#").split(delimiter)
        columns_of_interest = columns.strip().split(",") if columns else header_line[1:]
        logging.info("    Got header: %s", header_line)
        logging.info("    Got columns: %s", columns_of_interest)
        for col in columns_of_interest:
            if col not in header_line:
                logging.error("    column %s is not present in headers", col)
                raise ValueError(f"column {col} is not present in headers.")
        for n, line in enumerate(_in, start=1):
            line = line.strip()
            if line and line[0] != "#":
                line_d = dict(zip(header_line, line.split(delimiter)))
                for category in columns_of_interest:
                    features = line_d.get(category, "").strip()
                    if features:
                        cat_d.setdefault(category, set()).update(features.split(","))

        logging.info("    Parsed %s entries.", n)

        logging.info("Building code map and dumping category and feature encodings.")
        code_map = {}
        feature_offset = 0

        with gzip.open(db_path + ".code_map.json.gz", "wt") as _map_out:
            for category, features in sorted(cat_d.items()):
                code_map[category] = {
                    "key": len(code_map),
                    "features": {
                        feature: (i + feature_offset) for i, feature in enumerate(sorted(features))
                    }
                }
                feature_offset += len(features)

                if db_session is not None:
                    db_category = db.Category(id=code_map[category]["key"], name=category)
                    db_session.add(db_category)
                    db_session.commit()

                for feature, fid in code_map[category]["features"].items():
                    if db_session is not None:
                        db_feature = db.Feature(id=fid, name=feature, category=db_category)
                        db_session.add(db_feature)

                if db_session is not None:
                    db_session.commit()

            json.dump(code_map, _map_out)

        return code_map, n


def process_annotations(input_data, db_session, code_map, header=None, columns=None, delimiter="\t", is_gzipped=False):
    logging.info("Second pass: Encoding sequence annotations")

    with (gzip.open if is_gzipped else open)(input_data, "rt") as _in:
        if header:
            _ = [next(_in) for _ in range(header - 1)]
        header_line = next(_in).strip().strip("#").split(delimiter)
        columns_of_interest = columns.strip().split(",") if columns else header_line[1:]
        for _, line in enumerate(_in, start=1):
            line = line.strip()
            if line and line[0] != "#":
                line = line.split(delimiter)
                line_d = dict(zip(header_line, line))
                encoded = []
                for category in columns_of_interest:
                    features = line_d.get(category, "").strip()
                    if features:
                        enc_category = code_map[category]['key']
                        enc_features = sorted(
                            code_map[category]['features'][feature]
                            for feature in features.split(",")
                        )
                        encoded.append((enc_category, ",".join(map(str, enc_features))))
                encoded = ";".join(f"{cat}={features}" for cat, features in sorted(encoded))

                db_sequence = db.AnnotatedSequence(
                    seqid=line[0],
                    featureid=None,
                    strand=None,
                    annotation_str=encoded
                )
                db_session.add(db_sequence)

        db_session.commit()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("db_path", type=str)
    ap.add_argument("input_data", type=str)
    ap.add_argument("--initialise_db", action="store_true")
    ap.add_argument("--columns", type=str)
    ap.add_argument("--code_map", type=str)
    ap.add_argument("--nseqs", type=int)
    ap.add_argument("--extract_map_only", action="store_true")
    ap.add_argument("--header", type=int)
    ap.add_argument("--delimiter", type=str, default="\t")
    args = ap.parse_args()

    engine, db_session = get_database(args.db_path) if not args.extract_map_only else (None, None)

    if args.initialise_db and not args.extract_map_only:
        initialise_db(engine)

    gz_magic = b"\x1f\x8b\x08"
    # pylint: disable=R1732,W0511
    is_gzipped = open(args.input_data, "rb").read(3).startswith(gz_magic)

    if args.code_map:
        with gzip.open(args.code_map, "rt") as _map_in:
            code_map = json.load(_map_in)
    else:
        code_map, _ = gather_category_and_feature_data(
            args.input_data, args.db_path,
            db_session=db_session, columns=args.columns, header=args.header, delimiter=args.delimiter,
            is_gzipped=is_gzipped,
        )

    if args.extract_map_only:
        return

    process_annotations(
        args.input_data,
        db_session,
        code_map,
        columns=args.columns,
        header=args.header,
        delimiter="\t",
        is_gzipped=is_gzipped,
    )

    # https://www.sqlite.org/wal.html
    # https://stackoverflow.com/questions/10325683/can-i-read-and-write-to-a-sqlite-database-concurrently-from-multiple-connections
    # concurrent read-access from more than 3 processes seems to be an issue
    with sqlite3.connect(args.db_path) as conn:
        cur = conn.cursor()
        cur.execute('PRAGMA journal_mode=wal')
        cur.fetchall()


if __name__ == "__main__":
    main()

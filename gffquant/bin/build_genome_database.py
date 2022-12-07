# pylint: disable=R0914,C0103
# pylint: disable=duplicate-code

""" module docstring """

import argparse
import gzip
import json
import logging
import sqlite3

from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

from gffquant.db import initialise_db
from gffquant.db.models import db
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


def gather_category_and_feature_data(input_data, db_path, db_session=None, fields=None):

    cat_d = {}
    valid_fields = set(fields)
    logging.info("First pass: gathering category and feature information.")

    n = 0
    n_genes = 0
    with open(input_data, "rt") as _in:
        for n, (line) in enumerate(_in, start=1):
            if line[0] == "#":
                continue
            line = line.strip().split("\t")
            if line[2] != "gene":
                continue
            n_genes += 1
            attribs = dict(
                item.split("=")
                for item in line[8].strip(";").split(";")
            )
            for field in valid_fields:
                cat_d.setdefault(field, set()).update(set(attribs.get(field, "").split(",")) - {''})

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

    return code_map, n_genes


def process_annotations(input_data, db_session, code_map, nseqs):
    logging.info("Second pass: Encoding sequence annotations")

    n_genes = 0
    with open(input_data, "rt") as _in:
        for i, (line) in enumerate(_in, start=1):
            if line[0] == "#":
                continue

            line = line.strip().split("\t")
            if line[2].lower() != "gene":
                continue

            n_genes += 1
            
            if i % 10000 == 0:
                db_session.commit()

            if nseqs is not None:
                logging.info("Processed %s entries. (%s%%)", nseqs, round(n_genes / nseqs * 100, 3))
            else:
                logging.info("Processed %s entries.", str(n_genes))            

            attribs = dict(
                item.split("=")
                for item in line[8].strip(";").split(";")
            )

            encoded = []
            for category in code_map:
                enc_category = code_map[category]['key']
                features = set(attribs.get(category, "").split(",")) - {''}
                enc_features = sorted(code_map[category]['features'][feature] for feature in features)
                encoded.append((enc_category, ",".join(map(str, enc_features))))
            encoded = ";".join(f"{cat}={features}" for cat, features in sorted(encoded))

            strand = line[6] if line[6] != "." else None


            db_sequence = db.AnnotatedSequence(
                seqid=line[0],
                featureid=attribs.get("ID"),
                strand=int(strand == "+") if strand is not None else None,
                start=int(line[3]),
                end=int(line[4]),
                annotation_str=encoded
            )
            db_session.add(db_sequence)

        db_session.commit()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("db_path", type=str)
    ap.add_argument("input_data", type=str)
    ap.add_argument("--initialise_db", action="store_true")
    ap.add_argument("--code_map", type=str)
    ap.add_argument("--nseqs", type=int)
    ap.add_argument("--extract_map_only", action="store_true")
    ap.add_argument("--fields", type=str)
    
    args = ap.parse_args()

    engine, db_session = get_database(args.db_path) if not args.extract_map_only else (None, None)

    if args.initialise_db and not args.extract_map_only:
        initialise_db(engine)

    nseqs = args.nseqs
    if args.code_map:
        with gzip.open(args.code_map, "rt") as _map_in:
            code_map = json.load(_map_in)
    else:
        code_map, nseqs = gather_category_and_feature_data(args.input_data, args.db_path, db_session=db_session, fields=args.fields.split(",") if args.fields else [])

    if args.extract_map_only:
        return

    process_annotations(args.input_data, db_session, code_map, nseqs)

    # https://www.sqlite.org/wal.html
    # https://stackoverflow.com/questions/10325683/can-i-read-and-write-to-a-sqlite-database-concurrently-from-multiple-connections
    # concurrent read-access from more than 3 processes seems to be an issue
    with sqlite3.connect(args.db_path) as conn:
        cur = conn.cursor()
        cur.execute('PRAGMA journal_mode=wal')
        cur.fetchall()


if __name__ == "__main__":
    main()

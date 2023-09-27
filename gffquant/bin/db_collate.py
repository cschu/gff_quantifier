# pylint: disable=C0103,W1203,R0914,R0915

""" db_collate """

import argparse
import logging
import os
import pathlib
import sqlite3
import time

import pandas as pd

from sqlalchemy import create_engine, insert
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", "-i", type=str, default=".")
    ap.add_argument("--output_prefix", "-o", type=str, default="collated")
    ap.add_argument("--column", "-c", type=str)
    ap.add_argument("--db_path", type=str, default="observations.sqlite3")

    args = ap.parse_args()

    outdir = os.path.dirname(args.output_prefix)
    if outdir and outdir != ".":
        pathlib.Path(outdir).mkdir(exist_ok=True, parents=True)

    dirpath, _, files = next(os.walk(args.input_dir))

    engine, db_session = get_database(args.db_path)
    initialise_db(engine)

    features_d = {}
    samples_d = {}

    t000 = t00 = time.time()
    for i, f in enumerate(sorted(files), start=1):
        logging.info(f"Processing file {i}/{len(files)}: {f}")

        t0 = time.time()
        fname = os.path.basename(f).replace(".txt.gz", "")
        *sample, category = fname.split(".")
        sample = ".".join(sample)

        sample_id = samples_d.setdefault(sample, len(samples_d))
        # db_session.add(
        # 	db.Sample(id=sample_id, name=sample)
        # )

        counts = pd.read_csv(
            os.path.join(dirpath, f),
            delimiter="\t",
            index_col=0,
            usecols=(("feature", "gene")[category == "gene_counts"], args.column),
        )

        observations = []
        for _, row in counts.iterrows():
            if row.name != "unannotated":
                feature_name = row.name
                feature_id = features_d.setdefault(feature_name, len(features_d))

                observations.append(
                    {
                        "value": float(row[0]),
                        "category_id": 0,
                        "sample_id": sample_id,
                        "feature_id": feature_id,
                    }
                )

                # db_observation = db.Observation(
                # 	metric=args.column,
                # 	value=float(row[0]),
                # 	category_id=0,
                # 	sample_id=sample_id,
                # 	feature_id=feature_id,
                # )
                # db_session.add(db_observation)

        # db.Observation.__table__.insert().execute(observations)
        db_session.execute(insert(db.Observation), observations)

        db_session.commit()
        logging.info(f"Finished loading {f} in {time.time() - t0}s.")

    logging.info(f"Finished loading all files in {time.time() - t00}s.")

    logging.info(f"Adding {len(features_d)} features to database...")
    t0 = time.time()
    # for feature_name, feature_id in features_d.items():
    # 	# logging.info(f"    {feature_name} -> {feature_id}")
    # 	db_feature = db.Feature(id=feature_id, name=feature_name, category_id=0)
    # 	db_session.add(db_feature)

    db_session.execute(
        insert(db.Feature),
        [
            {
                "id": feature_id,
                "name": feature_name,
                "category_id": 0,
            }
            for feature_name, feature_id in features_d.items()
        ]
    )

    logging.info(f"Finished in {time.time() - t0}s.")

    logging.info(f"Adding {len(samples_d)} samples to database...")
    t0 = time.time()

    db_session.execute(
        insert(db.Sample),
        [
            {
                "id": sample_id,
                "name": sample_name,
            }
            for sample_name, sample_id in samples_d.items()
        ]
    )

    logging.info(f"Finished in {time.time() - t0}s.")

    db_session.commit()

    features_d.clear()

    logging.info("Converting database to count matrix...")
    t0 = time.time()
    con = sqlite3.connect(args.db_path)
    df = pd.read_sql_query(
        "select sample.name as sample_id, "
        "feature.name as feature_id, "
        "observation.value "
        "from observation "
        "join sample on sample.id = observation.sample_id "
        "join feature on feature.id = observation.feature_id",
        con
    )

    df = df.pivot(index="feature_id", columns="sample_id")
    df.columns = [x[1] for x in df.columns]
    df.index.name = "feature"
    logging.info(f"Finished in {time.time() - t0}s.")
    logging.info("Saving database...")
    t0 = time.time()
    df.to_csv(
        f"{args.output_prefix}.{category}.{args.column}.txt.gz",
        sep="\t",
        na_rep="NA",
    )
    logging.info(f"Finished in {time.time() - t0}s.")

    logging.info(f"Collation complete in {time.time() - t000}s.")


if __name__ == "__main__":
    main()

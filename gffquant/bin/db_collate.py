import argparse
import csv
import gzip
import os
import sqlite3

import pandas as pd

from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

from ..db import initialise_db
from ..db.models import db
from ..db.models.meta import Base


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

	dirpath, _, files = next(os.walk(args.input_dir))

	engine, db_session = get_database(args.db_path)
	initialise_db(engine)


	for f in files:
		fname = os.path.basename(f).replace(".txt.gz", "")
		*sample, category = fname.split(".")
		sample = ".".join(sample)

		db_category = db_session.query(db.Category)\
			.filter(db.Category.name == category).one_or_none()
		if db_category is None:
			db_category = db.Category(name=category)
			# category_id = db_category.id
			db_session.add(db_category)
			db_session.commit()
		# else:
		# 	category_id = db_category.id
		

		db_sample = db_session.query(db.Sample)\
			.filter(db.Sample.name == sample).one_or_none()
		if db_sample is None:
			db_sample = db.Sample(name=sample)
			# sample_id = db_sample.id
			db_session.add(db_sample)
			db_session.commit()
		# else:
		# 	sample_id = db_sample.id		
		
		f_open = gzip.open if f.endswith(".gz") else open

		with f_open(os.path.join(dirpath, f), "rt") as _in:
			for row in csv.DictReader(_in, delimiter="\t"):
				if not "unannotated" in row:
					db_feature = db_session.query(db.Feature)\
						.filter(db.Feature.name == row["feature"]).one_or_none()
					if db_feature is None:
						db_feature = db.Feature(name=row["feature"], category_id=db_category.id)
						feature_id = db_feature.id
						db_session.add(db_feature)
						db_session.commit()
					# else:
					# 	feature_id = db_feature.id
					db_observation = db.Observation(
						metric=args.column,
						value=float(row[args.column]),
						category_id=db_category.id,
						sample_id=db_sample.id,
						feature_id=db_feature.id,
					)
					db_session.add(db_observation)
					db_session.commit()


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
	df.index.name = ""
	df.to_csv(
		f"{args.output_prefix}.{category}.{args.column}.txt.gz",
		sep="\t",
		na_rep="NA",
	)
	# collated.CAZy.uniq_scaled.txt.gz


	
	
if __name__ == "__main__":
	main()
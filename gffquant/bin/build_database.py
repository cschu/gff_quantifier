import argparse
import csv
import gzip
import logging


logging.basicConfig(
    # filename=filename_log,
    # filemode='w',
    level=logging.INFO,
    format='[%(asctime)s] %(message)s'
)

#from sqlalchemy import create_engine
#from sqlalchemy.orm import scoped_session, sessionmaker



# from ..db.models.meta import Base
from gffquant.db import get_database, initialise_db

from gffquant.db.models import db

from gffquant.db.gff_dbm import GffDatabaseManager


# # def get_database(db_path):
# # 	engine = create_engine(f"sqlite:///{db_path}")

# # 	db_session = scoped_session(
# # 		sessionmaker(
# # 			autocommit=False,
# # 			autoflush=False,
# # 			enable_baked_queries=True,
# # 			bind=engine
# # 		)
# # 	)

# # 	Base.query = db_session.query_property()

# # 	return engine, db_session

# # def initialise_db(engine):
# 	Base.metadata.create_all(bind=engine)



def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("db_path", type=str)
	ap.add_argument("input_data", type=str)
	ap.add_argument("--initialise_db", action="store_true")
	ap.add_argument(
        "--emapper_version",
        type=str,
        default="v2",
        choices=("v1", "v2"),
    )
	args = ap.parse_args()

	engine, db_session = get_database(args.db_path)

	if args.initialise_db:
		initialise_db(engine)

	cat_d = {}

	logging.info("First pass: gathering category and feature information.")
	gffdbm = GffDatabaseManager(args.input_data, "genes", emapper_version=args.emapper_version)

	n = 0	
	for n, (ref, region_annotation) in enumerate(gffdbm.iterate(bufsize=4000000000), start=1):
		for category, features in region_annotation[1:]:
			cat_d.setdefault(category, set()).update(features)

	logging.info(f"    Parsed {n} entries.")
	
	logging.info("Building code map and dumping category and feature encodings.")
	code_map = {}
	feature_offset = 0

	with gzip.open(args.db_path + ".category.map.gz", "wt") as cat_out, gzip.open(args.db_path + ".feature.map.gz", "wt") as feat_out:
		
		for category, features in sorted(cat_d.items()):
			code_map[category] = {
				"key": len(code_map),
				"features": {
					feature: (i + feature_offset) for i, feature in enumerate(sorted(features))
				}
			}
			feature_offset += len(features)
	
			print(category, code_map[category]["key"], sep="\t", file=cat_out)
	
			db_category = db.Category(id=code_map[category]["key"], name=category)
			db_session.add(db_category)
			db_session.commit()
			
			for feature, fid in code_map[category]["features"].items():
	
				print(feature, fid, sep="\t", file=feat_out)
				
				db_feature = db.Feature(id=fid, name=feature, category=db_category)
				db_session.add(db_feature)
			
			db_session.commit()

	
	logging.info("Second pass: Encoding sequence annotations")
	gffdbm = GffDatabaseManager(args.input_data, "genes", emapper_version=args.emapper_version)
	for i, (ref, region_annotation) in enumerate(gffdbm.iterate(bufsize=4000000000), start=1):
		if i % 10000 == 0:
			db_session.commit()
			logging.info(f"Processed {i} entries. ({i/n * 100:.03f}%)")

		encoded = []
		for category, features in region_annotation[1:]:
			enc_category = code_map[category]['key']
			enc_features = sorted(code_map[category]['features'][feature] for feature in features)
			encoded.append((enc_category, ",".join(map(str, enc_features))))
		encoded = ";".join(f"{cat}={features}" for cat, features in sorted(encoded))

		strand, gene_id = region_annotation[0]
		db_sequence = db.AnnotatedSequence(
			seqid=ref if gene_id is None else gene_id,
			contig=ref if gene_id is not None else None,
			strand=int(strand == "+") if strand is not None else None,
			annotation_str=encoded
		)
		db_session.add(db_sequence)
	db_session.commit()


if __name__ == "__main__":
	main()

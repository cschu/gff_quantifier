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
	for n, (ref, region_annotation) in enumerate(gffdbm.iterate(), start=1):
		for category, features in region_annotation[1:]:
			cat_d.setdefault(category, set()).update(features)

	logging.info(f"    Parsed {n} entries.")
	
	logging.info("Building code map and dumping category and feature encodings.")
	code_map = {}
	feature_offset = 0
	for category, features in sorted(cat_d.items()):
		code_map[category] = {
			"key": len(code_map),
			"features": {
				feature: (i + feature_offset) for i, feature in enumerate(sorted(features))
			}
		}
		feature_offset += len(features)

		db_category = db.Category(id=code_map[category]["key"], name=category)
		db_session.add(db_category)
		db_session.commit()
		
		for feature, fid in code_map[category]["features"].items():
			db_feature = db.Feature(id=fid, name=feature, category=db_category)
			db_session.add(db_feature)
		
		db_session.commit()

	
	logging.info("Second pass: Encoding sequence annotations")
	gffdbm = GffDatabaseManager(args.input_data, "genes", emapper_version=args.emapper_version)
	for i, (ref, region_annotation) in enumerate(gffdbm.iterate(), start=1):
		if i % 10000 == 0:
			logging.info(f"Processed {i} entries. ({i/n * 100:.03f}%)")
		encoded = []
		for category, features in region_annotation[1:]:
			enc_category = code_map[category]['key']
			enc_features = sorted(code_map[category]['features'][feature] for feature in features)
			encoded.append((enc_category, ",".join(map(str, enc_features))))
		encoded = ";".join(f"{cat}={features}" for cat, features in sorted(encoded))

		# print(ref, encoded, sep="\t")

		strand, gene_id = region_annotation[0]
		db_sequence = db.AnnotatedSequence(
			seqid=ref if gene_id is None else gene_id,
			contig=ref if gene_id is not None else None,
			strand=int(strand == "+") if strand is not None else None,
			annotation_str=encoded
		)
		db_session.add(db_sequence)
		db_session.commit()

	return None
	for ref, region_annotation in gffdbm.iterate():
		# yield line[self.emapper_format.query_field], (
        #                        ("strand", None),
        #                    ) + tuple(categories)
		strand, gene_id = region_annotation[0]
		db_sequence = db.Sequence(
			seqid=ref if gene_id is None else gene_id,
			contig=ref if gene_id is not None else None,
			strand=int(strand == "+") if strand is not None else None
		)
		db_session.add(db_sequence)
		db_session.commit()

		for category, features in region_annotation[1:]:

			db_category = db_session.query(db.FunctionalCategory).filter(db.FunctionalCategory.name == category).one_or_none()
			if db_category is None:
				db_category = db.FunctionalCategory(name=category)
				db_session.add(db_category)
				db_session.commit()

			for feature in features:
				db_feature = db_session.query(db.FunctionalFeature).filter(db.FunctionalFeature.name == feature).one_or_none()
				if db_feature is None:
					db_feature = db.FunctionalFeature(name=feature, category_id=db_category.id)
					db_session.add(db_feature)
					db_session.commit()

				db_annotation = db.Annotation(feature_id=db_feature.id, sequence_id=db_sequence.id)
				db_session.add(db_annotation)
				db_session.commit()

			# print(category, features)


		# sequence = db.Sequence(seqid=row["query_name"], contig=row["query_name"])
		# db_session.add(sequence)
		# db_session.commit()

		"""
			#query_name	seed_eggNOG_ortholog	seed_ortholog_value	seed_ortholog_score	Predicted_taxonomic_group	Predicted_protein_name	Gene_Ontology_terms	EC_number	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TCCAZy	BiGG_Reaction	tax_scope	eggNOG_OGs	bestOG	COG_Functional_Category	eggNOG_free_text
		"""
			
		








if __name__ == "__main__":
	main()

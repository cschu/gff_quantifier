# flake8: noqa

""" database importers """

from .annstr_db_importer import AnnstrDatabaseImporter
from .custom_database_importer import GqCustomDatabaseImporter
from .database_importer import GqDatabaseImporter
from .gene_database_importer import GqGeneDatabaseImporter
from .small_genome_database_importer import SmallGenomeDatabaseImporter
from .small_database_importer import SmallDatabaseImporter


def extract_features(columns):
	annotation = [
		(category, tuple(features.split(",")))
		for category, features in columns.items()
		if features and features != "-" and category != "COG_category"
	]

	cog_category = columns.get("COG_category")
	if cog_category and cog_category != "-":
		cog_category = cog_category.replace(",", "")
		if len(cog_category) > 1:
			annotation.append(("COG_category_composite", cog_category))
		annotation.append(("COG_category", tuple(cog_category)))

	return annotation
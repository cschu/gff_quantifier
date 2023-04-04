# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys

from gffquant import __version__ as gffquant_version
here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

	name="gffquant"
	version = gffquant_version

	if sys.version_info.major != 3:
		raise EnvironmentError("""{toolname} is a python module that requires python3, and is not compatible with python2.""".format(toolname=name))

	setup(
		name=name,
		version=version,
		description=description,
		long_description=long_description,
		url="https://github.com/cschu/gff_quantifier",
		author="Christian Schudoma",
		author_email="christian.schudoma@embl.de",
		license="MIT",
		classifiers=[
			"Development Status :: 4 - Beta",
			"Topic :: Scientific Engineering :: Bio/Informatics",
			"License :: OSI Approved :: MIT License",
			"Operating System :: POSIX :: Linux",
			"Programming Language :: Python :: 3.7"
		],
		zip_safe=False,
		keywords="large reference data genomic feature quantification",
		packages=find_packages(exclude=["test"]),
		#scripts=["util/gff_indexer.py"],
		install_requires=[
			"intervaltree",
			"numpy",
			"pandas",
			"sqlalchemy",
			"pysam",
		],
		entry_points={
			"console_scripts": [
				"gffquant=gffquant.__main__:main",
				"gffindex=gffquant.gff_indexer:main",
				"collate_counts=gffquant.bin.collate_counts:main",
				"split_table=gffquant.bin.split_table:main",
				"build_gene_database=gffquant.bin.build_gene_database:main",
				"build_domain_database=gffquant.bin.build_domain_database:main",
				"build_bed_database=gffquant.bin.build_bed_database:main",
				"build_custom_database=gffquant.bin.build_custom_database:main",
				"collate_studies=gffquant.bin.collate_studies:main",
				"db_collate=gffquant.bin.db_collate:main",
			],
		},
		package_data={},
		include_package_data=True,
		data_files=[],
	)

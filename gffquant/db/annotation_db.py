# pylint: disable=C0103,R0903,W1514

""" module docstring """

from nis import cat
import sys
import gzip
from functools import lru_cache

from intervaltree import IntervalTree
from gffquant import annotation

from gffquant.db import get_database
from gffquant.db.models import db


"""
yield line[self.emapper_format.query_field], (
                                ("strand", None),
                            ) + tuple(categories)
if features:
	categories.append((category, features))
"""


class AnnotationDatabaseManager:
	def __init__(self, db_path):
		_, self.dbsession = get_database(db_path)

	def query_sequence(self, seqid):
		db_sequence = self.dbsession.query(db.AnnotatedSequence).filter(db.AnnotatedSequence.seqid == seqid).one_or_none()
		if db_sequence is None:
			return None
		categories = tuple()
		for cat_features in db_sequence.annotation_str.strip().split(";"):
			category, features = cat_features.split("=")
			categories += ((category, tuple(feature.strip() for feature in features.split(",") if feature.strip())),)
		return categories

	def query_feature(self, feature_id):
		return self.dbsession.query(db.Feature).filter(db.Feature.id == feature_id).join(db.Category, db.Feature.category_id == db.Category.id).one_or_none()
	
	def query_category(self, category_id):
		return self.dbsession.query(db.Category).filter(db.Category.id == category_id).one_or_none()

	def get_overlaps(self, ref, start, end):
		raise NotImplementedError
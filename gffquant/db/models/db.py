from ast import For
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, Date, ForeignKey, DateTime, Boolean, null, Table
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from .meta import Base


class AnnotatedSequence(Base):
	__tablename__ = "annotatedsequence"

	id = Column(Integer, primary_key=True)
	seqid = Column(String, index=True)
	contig = Column(String)
	start = Column(Integer)
	end = Column(Integer)
	strand = Column(Integer)

	annotation_str = Column(String)


class Category(Base):
	__tablename__ = "category"

	id = Column(Integer, primary_key=True)
	name = Column(String)

	features = relationship("Feature", back_populates="category")


class Feature(Base):
	__tablename__ = "feature"

	id = Column(Integer, primary_key=True)
	name = Column(String)

	category_id = Column(Integer, ForeignKey("category.id"))

	category = relationship("Category", back_populates="features")

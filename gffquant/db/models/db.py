from ast import For
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, Date, ForeignKey, DateTime, Boolean, null, Table
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from .meta import Base


annotation_feature_table = Table(
	"sequence_annotation",
	Base.metadata,
	Column("annotation_id", Integer, ForeignKey("annotation.id")),
	Column("feature_id", Integer, ForeignKey("functionalfeature.id"))

)


class Sequence(Base):
	__tablename__ = "sequence"

	id = Column(Integer, primary_key=True)
	seqid = Column(String)
	contig = Column(String)
	start = Column(Integer)
	end = Column(Integer)
	strand = Column(Integer)

	# relationships
	annotations = relationship("Annotation", back_populates="sequence")


class Annotation(Base):
	__tablename__ = "annotation"

	id = Column(Integer, primary_key=True)
	feature_id = Column(String, ForeignKey("functionalfeature.id"))

	sequence_id = Column(Integer, ForeignKey("sequence.id"))
	
	# relationships
	sequence = relationship("Sequence", back_populates="annotations")
	feature = relationship("FunctionalFeature", back_populates="annotations", secondary=annotation_feature_table)


class FunctionalCategory(Base):
	__tablename__ = "functionalcategory"

	id = Column(Integer, primary_key=True)
	name = Column(String)

	# relationships
	features = relationship("FunctionalFeature", back_populates="category")


class FunctionalFeature(Base):
	__tablename__ = "functionalfeature"

	id = Column(Integer, primary_key=True)
	name = Column(String)

	category_id = Column(Integer, ForeignKey("functionalcategory.id"))
	annotation_id = Column(Integer, ForeignKey("annotation.id"))

	# relationships
	category = relationship("FunctionalCategory", back_populates="features")
	annotations = relationship("Annotation", back_populates="feature", secondary=annotation_feature_table)
# pylint: disable=R0903

""" module docstring """

from sqlalchemy import Column, Integer, String, ForeignKey, Float
from sqlalchemy.orm import relationship

from .meta import Base


class AnnotatedSequence(Base):
    __tablename__ = "annotatedsequence"

    id = Column(Integer, primary_key=True)
    seqid = Column(String, index=True)
    featureid = Column(String)
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


class Sample(Base):
    __tablename__ = "sample"

    id = Column(Integer, primary_key=True)
    name = Column(String)

    observations = relationship("Observation", back_populates="sample")


class Observation(Base):
    __tablename__ = "observation"

    id = Column(Integer, primary_key=True)
    metric = Column(String)
    value = Column(Float)

    feature_id = Column(Integer, ForeignKey("feature.id"))
    category_id = Column(Integer, ForeignKey("category.id"))
    sample_id = Column(Integer, ForeignKey("sample.id"))

    sample = relationship("Sample", back_populates="observations")

# pylint: disable=R0903

""" module docstring """

from sqlalchemy import Column, Integer, String, ForeignKey
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

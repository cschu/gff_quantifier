"""module docstring"""

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.schema import MetaData

metadata = MetaData()
Base = declarative_base(metadata=metadata)

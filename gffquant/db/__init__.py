# pylint: disable=C0301
""" module docstring """

import sqlite3

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from .models import meta


def get_database(db_path, in_memory=True):

    if in_memory:
        source = sqlite3.connect(db_path)

        # https://stackoverflow.com/questions/68286690/copy-an-sqlite-database-into-memory-using-sqlalchemy-for-testing-flask-app !!!
        engine = create_engine("sqlite://", poolclass=StaticPool, connect_args={'check_same_thread': False},)
        conn = engine.raw_connection().connection
        source.backup(conn)

        meta.Base.metadata.create_all(engine)

    else:
        engine = create_engine(f"sqlite:///{db_path}")

    session = sessionmaker(bind=engine)
    db_session = session()

    return engine, db_session


def initialise_db(engine):
    meta.Base.metadata.create_all(bind=engine)

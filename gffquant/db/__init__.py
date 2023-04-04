# pylint: disable=C0301
""" module docstring """

import logging
import sqlite3

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import StaticPool

from .models import meta


logger = logging.getLogger(__name__)


def get_writable_database():
    # https://stackoverflow.com/questions/68286690/copy-an-sqlite-database-into-memory-using-sqlalchemy-for-testing-flask-app !!!
    engine = create_engine("sqlite://", poolclass=StaticPool, connect_args={'check_same_thread': False},)
    meta.Base.metadata.create_all(engine)

    session = sessionmaker(bind=engine)
    db_session = session()

    return engine, db_session


def get_database(db_path, in_memory=True):

    if in_memory:
        source = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)

        # https://stackoverflow.com/questions/68286690/copy-an-sqlite-database-into-memory-using-sqlalchemy-for-testing-flask-app !!!
        engine = create_engine("sqlite://", poolclass=StaticPool, connect_args={'check_same_thread': False},)
        conn = engine.raw_connection().connection
        logger.info("Loading database into memory...")
        source.backup(conn)
        logger.info("Finished loading database.")

        meta.Base.metadata.create_all(engine)

    else:
        engine = create_engine(f"sqlite:///{db_path}")

    session = sessionmaker(bind=engine)
    db_session = session()

    return engine, db_session


def initialise_db(engine):
    meta.Base.metadata.create_all(bind=engine)

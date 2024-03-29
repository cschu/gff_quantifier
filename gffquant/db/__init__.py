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


def initialise_db(engine):
    meta.Base.metadata.create_all(bind=engine)


def get_database(db_path, in_memory=True):

    if in_memory:
        source = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)

        # https://stackoverflow.com/questions/68286690/copy-an-sqlite-database-into-memory-using-sqlalchemy-for-testing-flask-app !!!
        engine = create_engine("sqlite://", poolclass=StaticPool, connect_args={'check_same_thread': False},)
        conn = engine.raw_connection().connection
        logger.info("Loading database into memory...")
        source.backup(conn)
        logger.info("Finished loading database.")

        # meta.Base.metadata.create_all(engine)
        initialise_db(engine)

    else:
        engine = create_engine(f"sqlite:///{db_path}")

    session = sessionmaker(bind=engine)
    db_session = session()

    return engine, db_session


# pylint: disable=C0103
def improve_concurrent_read_access(db):
    # https://www.sqlite.org/wal.html
    # https://stackoverflow.com/questions/10325683/can-i-read-and-write-to-a-sqlite-database-concurrently-from-multiple-connections
    # concurrent read-access from more than 3 processes seems to be an issue
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute('PRAGMA journal_mode=wal')
        cur.fetchall()

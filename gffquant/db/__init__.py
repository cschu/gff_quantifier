from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

from .models import meta



def get_database(db_path):
	engine = create_engine(f"sqlite:///{db_path}")

	db_session = scoped_session(
		sessionmaker(
			autocommit=False,
			autoflush=False,
			enable_baked_queries=True,
			bind=engine
		)
	)

	meta.Base.query = db_session.query_property()

	return engine, db_session

def initialise_db(engine):
	meta.Base.metadata.create_all(bind=engine)
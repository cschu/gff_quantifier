import sqlite3

from sqlalchemy import create_engine, text
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy import create_engine, MetaData, Table,Column,Integer,select
from sqlalchemy.orm import mapper, sessionmaker
from sqlite3 import dbapi2 as sqlite
from sqlalchemy.engine.reflection import Inspector
from sqlalchemy.pool import StaticPool

from .models import meta



def get_database(db_path, in_memory=True):

	if in_memory:
		source = sqlite3.connect(db_path)
		#conn = sqlite3.connect('file::memory:?cache=shared', uri=True)
		#source.backup(conn)
		#source.close()


		# # creating a random name for the temporary memory DB
		# sqlite_shared_name = "test_db_{}".format(
        #     random.sample(string.ascii_letters, k=4)
        # )

		# create_engine(
   		# 	"sqlite:///file:{}?mode=memory&cache=shared&uri=true".format(
        # 	sqlite_shared_name))


		engine = create_engine("sqlite://", poolclass=StaticPool, connect_args={'check_same_thread':False},)
		#engine.execute("ATTACH DATABASE 'file::memory:?cache=shared' AS 'my_db'")
		# meta.Base.metadata.drop_all(engine)

		conn = engine.raw_connection().connection
		source.backup(conn)


		meta.Base.metadata.create_all(engine)


		# metadata = MetaData(engine)


		#inspector = Inspector.from_engine(engine)
		#print(inspector.get_table_names())



		#moz_bookmarks = Table('table_a', metadata,Column("id", Integer, primary_key=True),schema='AA', autoload=True)
		#mapper(Bookmarks, moz_bookmarks)
		#moz_bookmarksB = Table('table_b', metadata,Column("id", Integer, primary_key=True),schema='BB', autoload=True)
		#mapper(BookmarksB, moz_bookmarksB)

		Session = sessionmaker(bind=engine)
		session = Session()

		#  meta.Base.query = session.query_property()

		return engine, session



		with engine.connect() as connection:
			connection.execute(text("ATTACH DATABASE ':memory:' AS my_db"))
			# connection.create_all()
		
		db_session = scoped_session(
			sessionmaker(
				autocommit=False,
				autoflush=False,
				enable_baked_queries=True,
				bind=engine
			)
		)

		meta.Base.query = db_session.query_property()

		#db_session.execute("ATTACH DATABASE ':memory:' AS my_db")
		#db_session.
				

	else:
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
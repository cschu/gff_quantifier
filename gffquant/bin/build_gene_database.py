# pylint: disable=R0914,C0103
# pylint: disable=duplicate-code

""" module docstring """

import argparse
import logging

from os.path import basename, splitext

from ..db import get_database, initialise_db, improve_concurrent_read_access
from ..db.importers import GqGeneDatabaseImporter

from .. import __tool__, __version__


logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s'
)

logger = logging.getLogger(__name__)


def main():
    ap = argparse.ArgumentParser(prog=f"{__tool__}:{splitext(basename(__file__))[0]}")
    ap.add_argument("db_path", type=str)
    ap.add_argument("input_data", type=str)
    ap.add_argument("--initialise_db", action="store_true")
    ap.add_argument("--code_map", type=str)
    ap.add_argument("--nseqs", type=int)
    ap.add_argument("--extract_map_only", action="store_true")
    ap.add_argument(
        "--emapper_version",
        type=str,
        default="v2",
        choices=("v1", "v2", "v2.1.2"),
    )
    ap.add_argument(
        "--version", "-v", action="version", version="%(prog)s " + __version__
    )
    args = ap.parse_args()

    engine, db_session = get_database(args.db_path, in_memory=False) if not args.extract_map_only else (None, None)

    if args.initialise_db and not args.extract_map_only:
        initialise_db(engine)

    _ = GqGeneDatabaseImporter(
        args.input_data,
        db_path=args.db_path,
        db_session=db_session,
        emapper_version=args.emapper_version,
    )

    improve_concurrent_read_access(args.db_path)


if __name__ == "__main__":
    main()

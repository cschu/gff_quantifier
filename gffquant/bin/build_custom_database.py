# pylint: disable=R0914,C0103,R0913,W1514
# pylint: disable=duplicate-code

""" module docstring """

import argparse
import logging

from os.path import basename, splitext

from ..db import get_database, initialise_db, improve_concurrent_read_access
from ..db.importers import AnnstrDatabaseImporter, GqCustomDatabaseImporter

from .. import __tool__, __version__


logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s'
)


def main():
    ap = argparse.ArgumentParser(prog=f"{__tool__}:{splitext(basename(__file__))[0]}")
    ap.add_argument("db_path", type=str)
    ap.add_argument("input_data", type=str)
    ap.add_argument("--initialise_db", action="store_true")
    ap.add_argument("--columns", type=str)
    ap.add_argument("--seq_column", type=str)
    ap.add_argument("--seqid_column", type=str)
    ap.add_argument("--code_map", type=str)
    ap.add_argument("--nseqs", type=int)
    ap.add_argument("--extract_map_only", action="store_true")
    ap.add_argument("--skip_header_lines", type=int, default=0)
    ap.add_argument("--header", type=str)
    ap.add_argument("--delimiter", type=str, default="\t")
    ap.add_argument("--dbtype", choices=("seq", "str"), default="seq")
    ap.add_argument(
        "--version", "-v", action="version", version="%(prog)s " + __version__
    )
    args = ap.parse_args()

    engine, db_session = get_database(args.db_path, in_memory=False) if not args.extract_map_only else (None, None)

    if args.initialise_db and not args.extract_map_only:
        initialise_db(engine)

    Importer = {"seq": GqCustomDatabaseImporter, "str": AnnstrDatabaseImporter}[args.dbtype]

    kwargs = {
        "seq_column": args.seq_column,
        "seqid_column": args.seqid_column,
    } if args.dbtype == "str" and args.seq_column else {}

    Importer(
        db_path=args.db_path,
        db_session=db_session,
        columns=args.columns.split(",") if args.columns is not None else None,
        skip_header_lines=args.skip_header_lines,
        header=args.header,
        delimiter=args.delimiter,
        **kwargs,
    ).build_database(
        args.input_data,
    )

    improve_concurrent_read_access(args.db_path)


if __name__ == "__main__":
    main()

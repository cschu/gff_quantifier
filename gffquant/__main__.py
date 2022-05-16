# pylint: disable=C0103,C0301

""" module docstring """

import logging
import os
import pathlib
import sys

from gffquant.feature_quantifier import FeatureQuantifier
from . import __version__
from .handle_args import handle_args


logger = logging.getLogger(__name__)


def main():

    args = handle_args(sys.argv[1:])

    print("Version:", __version__)
    print("Command:", os.path.basename(sys.argv[0]), *sys.argv[1:])

    if not os.path.exists(args.bam_file):
        raise ValueError("bam file does not exist", args.bam_file)
    if not os.path.exists(args.annotation_db):
        raise ValueError("annotation database does not exist", args.annotation_db)

    db_index = None
    if args.mode == "genome":
        db_index = args.annotation_db + ".index"
        if not os.path.exists(db_index):
            raise ValueError(
                f"gff index '{db_index}' does not exist "
                f"(please generate index with `gffindex {args.annotation_db}`)"
            )

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
        )

    fq = FeatureQuantifier(
        db=args.annotation_db,
        out_prefix=args.out_prefix,
        ambig_mode=args.ambig_mode,
        reference_type=args.mode,
        strand_specific=args.strand_specific,
        debugmode=args.debug,
    )

    fq.process_bamfile(
        args.bam_file, min_identity=args.min_identity, min_seqlen=args.min_seqlen
    )


if __name__ == "__main__":
    main()

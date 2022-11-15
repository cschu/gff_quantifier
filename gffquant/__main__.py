# pylint: disable=C0103,C0301

""" module docstring """

import logging
import os
import pathlib
import sys

# pylint: disable=W0611
from gffquant.profilers import GeneQuantifier, RegionQuantifier
from . import __version__
from .handle_args import handle_args


logger = logging.getLogger(__name__)


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
    logger.info("Command: %s %s", os.path.basename(sys.argv[0]), " ".join(sys.argv[1:]))

    if args.bam_file != "-" and not os.path.exists(args.bam_file):
        raise ValueError("bam file does not exist", args.bam_file)
    if not os.path.exists(args.annotation_db):
        raise ValueError("annotation database does not exist", args.annotation_db)

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
        )

    if args.mode in ("gene", "genes"):

        fq = GeneQuantifier(
            db=args.annotation_db,
            out_prefix=args.out_prefix,
            ambig_mode=args.ambig_mode,
            strand_specific=args.strand_specific,
            paired_end_count=args.paired_end_count,
            unmarked_orphans=args.unmarked_orphans,
        )

    else:

        fq = RegionQuantifier(
            db=args.annotation_db,
            out_prefix=args.out_prefix,
            ambig_mode=args.ambig_mode,
            strand_specific=args.strand_specific,
            calc_coverage=args.calc_coverage,
            paired_end_count=args.paired_end_count,
            unmarked_orphans=args.unmarked_orphans,
            reference_type=args.mode,
        )

    fq.process_bamfile(
        args.bam_file, aln_format=args.format, min_identity=args.min_identity, min_seqlen=args.min_seqlen
    )


if __name__ == "__main__":
    main()

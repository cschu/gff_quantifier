# pylint: disable=C0301,C0103
""" module docstring """

import argparse
import logging
import textwrap

from . import __version__


def handle_args(args):

    log_ap = argparse.ArgumentParser(prog="gffquant", add_help=False)
    log_ap.add_argument("-l", "--log_level", type=int, choices=range(1, 5), default=logging.INFO)
    log_args, _ = log_ap.parse_known_args(args)

    try:
        logging.basicConfig(
            level=log_args.log_level,
            format='[%(asctime)s] %(message)s'
        )
    except ValueError as invalid_loglevel_err:
        raise ValueError(f"Invalid log level: {log_args.log_level}") from invalid_loglevel_err

    ap = argparse.ArgumentParser(
        prog="gffquant",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=(log_ap,),
    )
    ap.add_argument(
        "annotation_db",
        type=str,
        help=textwrap.dedent(
            """\
			Path to a text file containing the reference annotation.
			The required type of file is determined by the --mode argument (gff3 or tsv)."""
        ),
    )
    ap.add_argument(
        "bam_file",
        type=str,
        help=textwrap.dedent(
            """\
			Path to a position-sorted bam file. Ambiguous alignments need to be flagged as secondary
			alignments with the same read id as their primary alignment.
			(e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0."""
        ),
    )
    ap.add_argument(
        "--mode",
        "-m",
        type=str,
        default="genome",
        choices=("genome", "genes", "gene", "domain"),
        help=textwrap.dedent(
            """\
			Run mode:"
			 - 'genome' counts reads aligned against contigs, which are annotated with a gff3 file.
				The gff3 needs to have been indexed with gffindex prior to the run.
			 - 'gene' counts reads aligned against gene sequences, which are annotated with a tab-separated file.
			 - 'genes' is an alias for the 'gene' mode
			 - 'domain' counts reads against domain annotations within gene sequences, which are annotated with a bed4 file."""
        ),
    )
    ap.add_argument(
        "--out_prefix",
        "-o",
        type=str,
        default="gffquant",
        help="Prefix for output files.",
    )
    ap.add_argument(
        "--ambig_mode",
        type=str,
        choices=("unique_only", "all1", "primary_only", "1overN"),
        default="unique_only",
        help=textwrap.dedent(
            """\
			Setting how ambiguous alignments should be treated. This setting mimics NGLess' behaviour.
			- 'unique_only' ignores any alignment flagged as ambiguous (MAPQ=0). This is the default setting.
			- 'all1' treats each alignment as unique (each ambiguous alignment contributes 1 count to features it aligns to.)
			- 'primary_only' takes the unique alignments and the primary and alignment of each ambiguous read group.
			- '1overN' each alignment contributes 1/(n=number of ambiguous alignments of the same read) counts to features it aligns to."""
        ),
    )

    ap.add_argument(
        "--strand_specific",
        action="store_true",
        help="Perform strand-specific counting for RNAseq reads. "
        "This flag is currently ignored for paired-end data.",
    )

    ap.add_argument(
        "--calc_coverage",
        action="store_true",
        help="Perform coverage calculations."
    )

    ap.add_argument(
        "--min_identity",
        type=float,
        default=0.97,
        help="Minimum sequence identity [n_match/length] for an alignment to be considered.",
    )

    ap.add_argument(
        "--min_seqlen",
        type=int,
        default=45,
        help="Minimum read length [bp] for an alignment to be considered.",
    )

    ap.add_argument(
        "--version", "-v", action="version", version="%(prog)s " + __version__
    )
    ap.add_argument("--debug", action="store_true")

    return ap.parse_args(args)

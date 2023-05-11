# pylint: disable=C0301,C0103
""" module docstring """

import argparse
import logging
import textwrap

from . import __version__
from . import __tool__


def handle_args(args):

    log_ap = argparse.ArgumentParser(prog=__tool__, add_help=False)
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
        prog=__tool__,
        formatter_class=argparse.RawTextHelpFormatter,
        parents=(log_ap,),
    )
    ap.add_argument(
        "annotation_db",
        type=str,
        help=textwrap.dedent(
            """\
            Path to an sqlite3 database containing the reference annotation.
			"""
        ),
    )
    ap.add_argument(
        "bam_file",
        type=str,
        help=textwrap.dedent(
            """\
            Path to a name-sorted sam or bam (s. --format) file. Ambiguous alignments need to be flagged as secondary
            alignments with the same read id as their primary alignment.
            (e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0.
            Input from STDOUT can be used with '-'."""
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
             - 'genome' counts reads aligned against contigs.
             - 'gene' (alias 'genes') counts reads aligned against gene sequences.
             - 'domain' counts reads against domain annotations within gene sequences."""
        ),
    )
    ap.add_argument(
        "--out_prefix",
        "-o",
        type=str,
        default=__tool__,
        help="Prefix for output files.",
    )
    ap.add_argument(
        "--ambig_mode",
        type=str,
        choices=("unique_only", "all1", "primary_only", "1overN"),
        default="unique_only",
        help=textwrap.dedent(
            """\
            Determines how ambiguous alignments should be treated. This setting mimics NGLess' behaviour.
            - 'unique_only' ignores any alignment flagged as ambiguous (MAPQ=0). This is the default setting.
            - 'all1' treats each alignment as unique (each ambiguous alignment contributes 1 count to features it aligns to.)
            - 'primary_only' takes the unique alignments and the primary alignment of each ambiguous read group.
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
        "--format",
        type=str,
        choices=("sam", "bam", "SAM", "BAM"),
        default="sam",
        help="Format of the alignment input. Supported: sam, bam.",
    )

    ap.add_argument(
        "--paired_end_count",
        type=int,
        choices=(1, 2),
        default=1,
        help="Paired-end count contribution: 0.5 / mate (1) or 1 / mate (2) [1]",
    )

    # orphan reads will not have flag 0x1 set
    ap.add_argument(
        "--unmarked_orphans",
        action="store_true",
        help="Ensure that alignments from unmarked orphan reads (from preprocessing) are properly accounted for.",
    )

    ap.add_argument(
        "--import_readcounts",
        type=int,
        help="Import externally derived readcounts to allow readcount-based normalisation for prefiltered bam files.",
    )

    ap.add_argument(
        "--restrict_metrics",
        type=str,
        help="Restrict reported metrics. Comma-separated list of `raw`, `lnorm`, `scaled`, `rpkm`.",
        default="raw,lnorm,scaled",
    )

    ap.add_argument(
        "--version", "-v", action="version", version="%(prog)s " + __version__
    )
    ap.add_argument("--debug", action="store_true")

    return ap.parse_args(args)

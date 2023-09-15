# pylint: disable=C0301,C0103,W1203
""" module docstring """

import argparse
import logging
import os
import textwrap

from . import __version__
from . import __tool__

from .ui.validation import check_bwa_index, check_minimap2_index


logger = logging.getLogger(__name__)


def validate_args(args):

    logger.info(f"args: {args.__dict__}")

    if not os.path.isfile(args.annotation_db):
        raise ValueError(f"Cannot find annotation db at `{args.annotation_db}`.")
    if (args.aligner == "bwa" and not check_bwa_index(args.reference)) or (args.aligner == "minimap" and not check_minimap2_index(args.reference)):
        raise ValueError(f"Cannot find reference index at `{args.reference}`.")

    has_fastq = any(
        map(
            lambda x: x is not None,
            (
                getattr(args, arg) for arg in ("reads1", "reads2", "singles", "orphans") if hasattr(args, arg)
            )
        )
    )

    if tuple(map(bool, (has_fastq, args.bam, args.sam))).count(True) != 1:
        raise ValueError(f"Need exactly one type of input: bam={bool(args.bam)} sam={bool(args.sam)} fastq={bool(has_fastq)}.")

    args.input_type = "fastq" if has_fastq else ("bam" if args.bam else "sam")

    if (args.reference or args.aligner) and not has_fastq:
        raise ValueError("--reference/--aligner are not needed with alignment input (bam, sam).")
    if bool(args.reference and args.aligner) != has_fastq:
        raise ValueError("--fastq requires --reference and --aligner to be set.")

    if args.restrict_metrics:
        restrict_metrics = set(args.restrict_metrics.split(","))
        invalid = restrict_metrics.difference(('raw', 'lnorm', 'scaled', 'rpkm'))
        if invalid:
            raise ValueError(f"Invalid column(s) in `--restrict_metrics`: {str(invalid)}")
        args.restrict_metrics = tuple(restrict_metrics)

    if os.path.isdir(os.path.dirname(args.out_prefix)) and not args.force_overwrite:
        raise ValueError(f"Output directory exists {os.path.dirname(args.out_prefix)}. Specify -f to overwrite.")

    return args


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
        "--db",
        dest="annotation_db",
        required=True,
        type=str,
        help=textwrap.dedent(
            """\
            Path to the reference annotation database.
			"""
        ),
    )

    ap.add_argument(
        "--db_coordinates",
        type=str,
        default="bed",
        choices=("bed", "hmmer"),
        help="Coordinate format for text-based annotation databases. bed=[start, end), hmmer=[start, end]"
    )

    ap.add_argument(
        "--db_separator",
        type=str,
        default="\t",
        help="Separator-character for text-based annotation databases."
    )

    # ap.add_argument(
    #     "--bed4",
    #     dest="annotation_bed4",
    #     type=str,
    #     help=textwrap.dedent(
    #         """\
    #         Path to a bed4 file containing the reference annotation.
    # 		"""
    #     ),
    # )

    ap.add_argument(
        "--reference",
        type=str,
        help=textwrap.dedent(
            """\
            Path to a (BWA, minimap) reference index.
			"""
        ),
    )

    ap.add_argument(
        "--aligner",
        type=str,
        help=textwrap.dedent(
            """\
            Select aligner to map fastq files against a reference index.
			"""
        ),
        choices=("bwa", "minimap2"),
    )

    ap.add_argument(
        "--keep_alignment_file",
        type=str,
        help="Save alignments in sam format to file."
    )
    # ap.add_argument(
    #     "input_files",
    #     type=str,
    #     nargs="+",
    #     help=textwrap.dedent(
    #         """\
    #         Path to a name-sorted sam or bam (s. --format) file. Ambiguous alignments need to be flagged as secondary
    #         alignments with the same read id as their primary alignment.
    #         (e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0.
    #         Input from STDOUT can be used with '-'."""
    #     ),
    # )
    ap.add_argument(
        "--bam",
        type=str,
        # nargs="+",
        help=textwrap.dedent(
            """\
            Path to a name-sorted BAM file. Ambiguous alignments need to be flagged as secondary
            alignments with the same read id as their primary alignment.
            (e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0.
            Input from STDIN can be specified with '-'."""
        ),
    )
    ap.add_argument(
        "--sam",
        type=str,
        # nargs="+",
        help=textwrap.dedent(
            """\
            Path to a name-sorted SAM file. Ambiguous alignments need to be flagged as secondary
            alignments with the same read id as their primary alignment.
            (e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0.
            Input from STDIN can be specified with '-'."""
        ),
    )

    ap.add_argument(
        "--fastq-r1",
        dest="reads1",
        nargs="*",
        type=str,
        help="A comma-delimited string of forward/R1 read fastq files."
    )

    ap.add_argument(
        "--fastq-r2",
        dest="reads2",
        nargs="*",
        type=str,
        help="A comma-delimited string of reverse/R2 read fastq files."
    )

    ap.add_argument(
        "--fastq-singles", "-s",
        dest="singles",
        nargs="*",
        type=str,
        help="A comma-delimited string of single-end read fastq files."
    )

    ap.add_argument(
        "--fastq-orphans",
        dest="orphans",
        nargs="*",
        type=str,
        help="A comma-delimited string of orphan read fastq files."
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
        default="1overN",
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
        "--force_overwrite",
        "-f",
        action="store_true",
        help="Overwrite existing output."
    )

    ap.add_argument(
        "--cpus_for_alignment", "-t",
        type=int, default=1,
        help="",
    )

    ap.add_argument(
        "--version", "-v", action="version", version="%(prog)s " + __version__
    )
    ap.add_argument("--debug", action="store_true")

    return validate_args(ap.parse_args(args))

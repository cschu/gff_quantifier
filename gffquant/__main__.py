import sys
import argparse
import os
import pathlib
import textwrap

from gffquant.feature_quantifier import FeatureQuantifier
from . import __version__

def main():
    ap = argparse.ArgumentParser(prog="gffquant", formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument(
        "annotation_db", type=str,
        help=textwrap.dedent("""\
            Path to a text file containing the reference annotation.
            The required type of file is determined by the --mode argument (gff3 or tsv)."""
        )
    )
    ap.add_argument(
        "bam_file", type=str,
        help=textwrap.dedent("""\
            Path to a position-sorted bam file. Ambiguous alignments need to be flagged as secondary
            alignments with the same read id as their primary alignment.
            (e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0."""
        )
    )
    ap.add_argument(
        "--mode", "-m", type=str, default="genome", choices=("genome", "genes", "gene", "domain"),
        help=textwrap.dedent("""\
            Run mode:"
             - 'genome' counts reads aligned against contigs, which are annotated with a gff3 file.
                The gff3 needs to have been indexed with gffindex prior to the run.
             - 'gene' counts reads aligned against gene sequences, which are annotated with a tab-separated file.
             - 'genes' is an alias for the 'gene' mode
             - 'domain' counts reads against domain annotations within gene sequences, which are annotated with a bed4 file."""
        )
    )
    ap.add_argument(
        "--out_prefix", "-o", type=str, default="gffquant",
        help="Prefix for output files."
    )
    ap.add_argument(
        "--ambig_mode", type=str, choices=("unique_only", "all1", "primary_only", "1overN"), default="unique_only",
        help=textwrap.dedent("""\
            Setting how ambiguous alignments should be treated. This setting mimics NGLess' behaviour.
            - 'unique_only' ignores any alignment flagged as ambiguous (MAPQ=0). This is the default setting.
            - 'all1' treats each alignment as unique (each ambiguous alignment contributes 1 count to features it aligns to.)
            - 'primary_only' takes the unique alignments and the primary and alignment of each ambiguous read group.
            - '1overN' each alignment contributes 1/(n=number of ambiguous alignments of the same read) counts to features it aligns to."""
        )
    )
    ap.add_argument(
        "--strand_specific", action="store_true",
        help="Perform strand-specific counting for RNAseq reads. This currently only works for single-end data. This flag is ignored for paired-end data."
    )
    ap.add_argument("--version", "-v", action="version", version="%(prog)s " + __version__)


    args = ap.parse_args()

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
            raise ValueError("gff index '{}' does not exist (please generate index with 'gffindex {}')".format(db_index, args.annotation_db))

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(exist_ok=True, parents=True)

    fq = FeatureQuantifier(
        args.annotation_db,
        db_index,
        out_prefix=args.out_prefix,
        ambig_mode=args.ambig_mode,
        do_overlap_detection=args.mode in ("genome", "domain"),
        strand_specific=args.strand_specific
    )

    fq.process_data(args.bam_file)


if __name__ == "__main__":
    main()

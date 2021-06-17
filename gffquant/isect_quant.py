import sys
import argparse
import os
import pathlib
import textwrap

from gffquant.feature_quantifier import FeatureQuantifier
from . import __version__


def main():
	ap = argparse.ArgumentParser(prog="isect_quant", formatter_class=argparse.RawTextHelpFormatter)
	ap.add_argument(
		"bed_file", type=str,
		help=textwrap.dedent("""\
			Path to a name-sorted bed file (bedtools intersect output).
			Ambiguous alignments need to be flagged as secondary
			alignments with the same read id as their primary alignment.
			(e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0."""
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

	if not os.path.exists(args.bed_file):
		raise ValueError("bed file does not exist", args.bed_file)

	if os.path.dirname(args.out_prefix):
		pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(exist_ok=True, parents=True)

	fq = FeatureQuantifier(
		out_prefix=args.out_prefix,
		ambig_mode=args.ambig_mode,
		do_overlap_detection=False,
		strand_specific=args.strand_specific
	)

	fq.process_bedfile(args.bed_file)


if __name__ == "__main__":
	main()

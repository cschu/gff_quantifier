import argparse
import os
import pathlib

from gffquant.feature_quantifier import FeatureQuantifier

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("gff_file", type=str, help="Path to a gff3 file containing the reference annotation. The gff3 needs to have been indexed with gffindex prior to the run.")
	ap.add_argument("bam_file", type=str,
		help="Path to a position-sorted bam file. Ambiguous alignments need to be flagged as secondary "
			"alignments with the same read id as their primary alignment."
			"(e.g. output from BWA mem -a). All alignments of an ambiguous group need to have MAPQ=0.")
	ap.add_argument("--out_prefix", "-o", type=str, default="gffquant", help="Prefix for output files.")
	ap.add_argument("--ambig_mode", type=str, choices=("unique_only", "all1", "1overN"), default="unique_only",
		help="Setting how ambiguous alignments should be treated. This setting mimics ngLess' behaviour.\n"
			"- 'unique_only' ignores any alignment flagged as ambiguous (MAPQ=0). This is the default setting.\n"
			"- 'all1' treats each alignment as unique (each ambiguous alignment contributes 1 count to features it aligns to)\n"
			"- '1overN' each alignment contributes 1/(n=number of ambiguous alignments of the same read) counts to features it aligns to."
	)
	ap.add_argument("--strand_specific", action="store_true",
		help="Perform strand-specific counting for RNAseq reads. This currently only works for single-end data. This flag is ignored for paired-end data.")
	args = ap.parse_args()

	gff_index = args.gff_file + ".index"

	if os.path.dirname(args.out_prefix):
		pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(exist_ok=True, parents=True)

	if not os.path.exists(args.gff_file):
		raise ValueError("gff database does not exist", args.gff_file)
	if not os.path.exists(gff_index):
		raise ValueError("gff index '{}' does not exist (please generate index with 'gffindex {}')".format(gff_index, args.gff_file))
	if not os.path.exists(args.bam_file):
		raise ValueError("bam file does not exist", args.bam_file)

	fq = FeatureQuantifier(
		args.gff_file,
		gff_index,
		out_prefix=args.out_prefix,
		ambig_mode=args.ambig_mode
	)

	fq.process_data(args.bam_file, strand_specific=args.strand_specific)


if __name__ == "__main__":
	main()

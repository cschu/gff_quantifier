import argparse
import os

from gffquant.feature_quantifier import FeatureQuantifier

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("gff_file", type=str)
	ap.add_argument("bam_file", type=str)
	ap.add_argument("out_prefix", type=str)
	ap.add_argument("--count_config", type=str)
	args = ap.parse_args()

	gff_index = args.gff_file + ".index"
	if not os.path.exists(args.gff_file):
		raise ValueError("gff database does not exist", args.gff_file)
	if not os.path.exists(gff_index):
		raise ValueError("gff index does not exist", gff_index)
	if not os.path.exists(args.bam_file):
		raise ValueError("bam file does not exist", args.bam_file)

	fq = FeatureQuantifier(
		args.gff_file,
		gff_index,
		count_config=args.count_config
	)
	fq.process_bam(args.bam_file, args.out_prefix)


if __name__ == "__main__":
	main()

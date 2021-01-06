import csv
import glob
import sys
import os
import argparse

def read_counts(count_file, columns=("uniq_scaled",)):
	counts = dict()
	with open(count_file) as count_stream:
		header = next(count_stream).strip().split("\t")
		counts["unannotated"] = next(count_stream).strip().split("\t")[1]
		for line in csv.DictReader(count_stream, fieldnames=header, delimiter="\t"):
			if line["subfeature"].startswith("#"):
				sf_counts = counts.setdefault(line["subfeature"][1:], dict())
			else:
				sf_counts[line["subfeature"]] = [line.get(col) for col in columns]
	return counts

def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("count_dir", type=str)
	ap.add_argument("--columns", type=str, default="uniq_scaled")
	args = ap.parse_args()

	columns = args.columns.split(",")
	counts = dict()

	features = dict()

	for i, f in enumerate(sorted(glob.glob(os.path.join(args.count_dir, "*.feature_counts.txt")))):
		key = os.path.basename(f).split(".")[0]
		counts[key] = read_counts(f, columns=columns)
		for feat, subfeat_d in counts[key].items():
			if feat != "unannotated":
				features.setdefault(feat, set()).update(subfeat_d.keys())
		if False: #i > 10:
			break

	for feature, subfeatures in features.items():
		with open(feature + ".summary.txt", "w") as _out:
			print("subfeature", *counts.keys(), sep="\t", file=_out)
			print("unannotated", *(fcounts["unannotated"] for _, fcounts in counts.items()), sep="\t", file=_out)
			for subfeature in sorted(subfeatures):
				subfeature_counts = [feature_counts.get(feature, dict()).get(subfeature, [0.0])[0] for sample, feature_counts in sorted(counts.items())]
				print(subfeature, *subfeature_counts, sep="\t", file=_out)
		
if __name__ == "__main__":
	main()		

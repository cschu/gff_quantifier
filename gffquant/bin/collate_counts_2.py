import argparse
import pathlib
import os

import pandas as pd


class FeatureCountCollator:
    def __init__(self, count_dir, prefix, column):
        self.count_dir = count_dir
        self.prefix = prefix
        self.column = column
        self.categories = {}
        self._collect_count_files()

    @staticmethod
    def is_valid_file(f):
        return all((
            f.endswith(".txt"),
            not f.endswith(".seqname.dist1.txt"),
            not f.endswith(".seqname.uniq.txt"),
            not f.endswith(".gene_counts.txt"),
            not f.endswith(".ambig_tmp.txt"),
        ))

    def _collect_count_files(self):
        all_files = []
        for pwd, _, files in os.walk(self.count_dir):
            all_files += (os.path.join(pwd, f) for f in files if FeatureCountCollator.is_valid_file(f))

        for f in all_files:
            sample, category = os.path.splitext(os.path.basename(f).replace(".txt", ""))
            self.categories.setdefault(category[1:], []).append((sample, f))

    def collate(self):
        for category, files in self.categories.items():
            print(category, self.column)
            self._collate_category(category, sorted(files))

    def _collate_category(self, category, files):
        with open(f"{self.prefix}.{category}.{self.column}.txt", "wt") as table_out:
            merged_tab = None
            for sample, fn in files:
                src_tab = pd.read_csv(fn, sep="\t", index_col=0)
                if merged_tab is None:
                    merged_tab = pd.DataFrame(index=src_tab.index)
                merged_tab = merged_tab.merge(src_tab[self.column], left_index=True, right_index=True, how="outer").rename(columns={self.column: sample})
                merged_tab[sample]["unannotated"] = src_tab["uniq_raw"]["unannotated"]
            merged_tab.to_csv(table_out, sep="\t", na_rep="NA")


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("count_dir", type=str)
	ap.add_argument("--out_prefix", "-o", type=str, default="./collated")
	args = ap.parse_args()

	outdir = os.path.dirname(args.out_prefix)
	if outdir and outdir != ".":
		pathlib.Path(outdir).mkdir(exist_ok=True, parents=True)

	FeatureCountCollator(args.count_dir, args.out_prefix, "uniq_scaled").collate()
	FeatureCountCollator(args.count_dir, args.out_prefix, "combined_scaled").collate()


if __name__ == "__main__":
	main()
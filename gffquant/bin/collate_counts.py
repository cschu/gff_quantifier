# pylint: disable=R0914,C0103,C0301,R0913
""" module docstring """

import argparse
import gzip
import os
import pathlib

import pandas as pd


class FeatureCountCollator:
    def __init__(self, count_dir, prefix, column, recursive=False, suffix=".txt.gz"):
        self.count_dir = count_dir
        self.prefix = prefix
        self.suffix = suffix
        self.column = column
        self.categories = {}
        self._collect_count_files(recursive=recursive)

    @staticmethod
    def is_valid_file(f, suffix):
        return all((
            f.endswith(suffix),
            not f.endswith(f".seqname.dist1{suffix}"),
            not f.endswith(f".seqname.uniq{suffix}"),
            not f.endswith(f".gene_counts{suffix}"),
            not f.endswith(f".ambig_tmp{suffix}"),
        ))

    def _collect_count_files(self, recursive=False):
        all_files = []
        for pwd, _, files in os.walk(self.count_dir):
            all_files += (os.path.join(pwd, f) for f in files if FeatureCountCollator.is_valid_file(f, self.suffix))
            if not recursive:
                break

        for f in all_files:
            sample, category = os.path.splitext(os.path.basename(f).replace(self.suffix, ""))
            self.categories.setdefault(category[1:], []).append((sample, f))

    def collate(self):
        for category, files in self.categories.items():
            print(category, self.column)
            self._collate_category(category, sorted(files))

    def _collate_category(self, category, files):
        table_file = f"{self.prefix}.{category}.{self.column}.txt.gz"
        index = set()
        for _, fn in files:
            with gzip.open(fn, "rt") as _in:
                index.update(row.strip().split("\t")[0] for row in _in if row.strip())
        merged_tab = pd.DataFrame(index=['unannotated'] + sorted(index.difference({'feature', 'unannotated'})))
        for sample, fn in files:
            src_tab = pd.read_csv(fn, sep="\t", index_col=0)
            merged_tab = merged_tab.merge(src_tab[self.column], left_index=True, right_index=True, how="outer")
            merged_tab.rename(columns={self.column: sample}, inplace=True)
            merged_tab[sample]["unannotated"] = src_tab["uniq_raw"]["unannotated"]
        merged_tab.to_csv(table_file, sep="\t", na_rep="NA", index_label="feature")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("count_dir", type=str)
    ap.add_argument("--out_prefix", "-o", type=str, default="./collated")
    ap.add_argument("--recursive", "-r", action="store_true")
    ap.add_argument("--column", "-c", type=str, choices=("uniq_raw", "uniq_lnorm", "uniq_scaled", "combined_raw", "combined_lnorm", "combined_scaled"), default="uniq_raw")
    args = ap.parse_args()

    outdir = os.path.dirname(args.out_prefix)
    if outdir and outdir != ".":
        pathlib.Path(outdir).mkdir(exist_ok=True, parents=True)

    FeatureCountCollator(args.count_dir, args.out_prefix, args.column, recursive=args.recursive).collate()


if __name__ == "__main__":
    main()

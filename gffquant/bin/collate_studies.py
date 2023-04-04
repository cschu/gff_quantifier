# pylint: disable=R0914,C0103,C0301,R0913,R0903
""" module docstring """

import argparse
import glob
import gzip
import os
import pathlib

import pandas as pd


class CollationCollator:
    def __init__(self, count_dir, prefix, pattern="*/collated.*.txt.gz"):
        self.files = {}
        self.prefix = prefix
        self.count_dir = count_dir
        self.pattern = pattern
        self._collect_count_files()

    def _collect_count_files(self):
        for fn in glob.glob(os.path.join(self.count_dir, self.pattern)):
            self.files.setdefault(os.path.basename(fn), []).append(fn)

    def collate(self):
        for fn, files in sorted(self.files.items()):
            basename = os.path.basename(self.prefix)
            if basename and self.prefix[-1] != ".":
                self.prefix += "."

            table_file = f"{self.prefix}{fn}"
            index = set()
            for fn in files:
                with gzip.open(fn, "rt") as _in:
                    index.update(row.strip().split("\t")[0] for row in _in if row.strip())
            merged_tab = pd.DataFrame(index=['unannotated'] + sorted(index.difference({'feature', 'unannotated'})))
            for fn in sorted(files):
                src_tab = pd.read_csv(fn, sep="\t", index_col=0)
                merged_tab = merged_tab.merge(src_tab, left_index=True, right_index=True, how="outer")
            merged_tab.to_csv(table_file, sep="\t", na_rep="NA", index_label="feature")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("count_dir", type=str)
    ap.add_argument("--out_prefix", "-o", type=str, default="./collated")
    ap.add_argument("--pattern", type=str, default="*/collated.*.txt.gz")
    args = ap.parse_args()

    outdir = os.path.dirname(args.out_prefix)
    if outdir and outdir != ".":
        pathlib.Path(outdir).mkdir(exist_ok=True, parents=True)

    CollationCollator(args.count_dir, args.out_prefix, pattern=args.pattern).collate()


if __name__ == "__main__":
    main()

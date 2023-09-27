# pylint: disable=R0914,C0103,C0301,R0913
""" module docstring """

import argparse
import gzip
import os
import pathlib
import time

import pandas as pd


def get_lines_from_chunks(f, bufsize=400000000):
    with (gzip.open if f.endswith(".gz") else open)(f, "r") as _in:
        tail = ""
        while True:
            chunk = "".join((tail, _in.read(bufsize).decode()))
            if not chunk:
                break
            chunk = chunk.split("\n")
            *chunk, tail = chunk
            for line in chunk:
                yield line
        if tail:
            yield tail


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
        def is_valid_suffix(f, suffix, wanted):
            has_suffix = f.endswith(suffix)
            return (has_suffix and wanted) or (not has_suffix and not wanted)

        return is_valid_suffix(f, suffix, True) and all(
            is_valid_suffix(f, f"{infix}{suffix}", False)
            for infix in (".seqname.dist1", ".seqname.uniq", ".ambig_tmp", "Counter", "domain.coverage")
        )

    def _collect_count_files(self, recursive=False):
        all_files = []
        for pwd, _, files in os.walk(self.count_dir):
            all_files += (os.path.join(pwd, f) for f in files if FeatureCountCollator.is_valid_file(f, self.suffix))
            if not recursive:
                break

        for f in all_files:
            sample, category = os.path.splitext(os.path.basename(f).replace(self.suffix, ""))
            self.categories.setdefault(category[1:], []).append((sample, f))

    def collate(self, index_file=None):
        for category, files in self.categories.items():
            print(category, self.column)
            if category == "aln_stats":
                self._collate_aln_stats(sorted(files))
            else:
                self._collate_category(category, sorted(files), index_file=index_file)

    # pylint: disable=W0613
    def _extract_index(self, files, feature_label):
        t00 = time.time()
        index = set()
        for _, fn in files:
            t0 = time.time()
            index.update(row.strip().split("\t")[0] for row in get_lines_from_chunks(fn) if row.strip())
            # index.update(pd.read_csv(fn, sep="\t", index_col=0, usecols=(feature_label,)).index)
            print(f"Updated index from {fn} in {time.time() - t0:.03f}s (|index|={len(index)}).", flush=True)
        print(f"Extracted feature column in {time.time() - t00:.03f}s.", flush=True)
        return index

    # pylint: disable=R0912,R0915
    def _collate_category(self, category, files, index_file=None):
        table_file = f"{self.prefix}.{category}.{self.column}.txt.gz"
        feature_label = ("feature", "gene")[category == "gene_counts"]
        print(index_file, index_file is None, os.path.isfile(str(index_file)))
        if index_file is None or not os.path.isfile(index_file):
            index = self._extract_index(files, feature_label)
            if index_file is not None:
                t0 = time.time()
                with gzip.open(index_file, "wt") as _out:
                    print(*index, sep="\n", file=_out)
                print(f"Saved index to `{index_file}` in {time.time() - t0:.03f}s.")
        else:
            print(f"Using existing index file `{index_file}`.", flush=True)
            index = {line.strip() for line in get_lines_from_chunks(index_file)}

        index_header = []
        if "total_reads" in index:
            index_header.append("total_reads")
        if "filtered_reads" in index:
            index_header.append("filtered_reads")
        if "category" in index:
            index_header.append("category")
        if "unannotated" in index:
            index_header.append("unannotated")

        index = index_header + sorted(
            index.difference(("unannotated", feature_label, "total_reads", "filtered_reads", "category"))
        )

        merged_tab = None

        t00 = time.time()
        col_buffer_size = 100
        col_buffer = []
        sum_loading_time, sum_time = 0.0, 0.0
        for i, (sample, fn) in enumerate(files):

            t0 = time.time()
            if i % col_buffer_size == 0:
                if i > 0:
                    col_buffer = [c for c in col_buffer if c is not None]
                    if col_buffer:
                        merged_tab = pd.concat(
                            ([merged_tab] if merged_tab is not None else []) + col_buffer,
                            axis=1
                        )

                        elapsed = time.time() - t0
                        sum_time += elapsed
                        print(f"Emptied column buffer in {elapsed:.03f}s.", flush=True)

                col_buffer = [None for _ in range(col_buffer_size)]

            t0 = time.time()
            try:
                data = pd.read_csv(
                    fn, sep="\t", index_col=0, usecols=(feature_label, self.column)
                )
            except (KeyError, ValueError):
                prev_col = self.column
                self.column = self.column.replace("combined_", "uniq_")
                try:
                    data = pd.read_csv(
                        fn, sep="\t", index_col=0, usecols=(feature_label, self.column)
                    )
                except (KeyError, ValueError) as err:
                    raise KeyError(f"Cannot find fallback-column `{self.column}` (from `{prev_col}` in file `{fn}`. Aborting.") from err

            col_buffer[i % col_buffer_size] = pd.Series(
                index=index,
                data=data[self.column],
                name=sample,
            )

            elapsed = time.time() - t0
            sum_time += elapsed
            sum_loading_time += elapsed

            print(
                f"Processed sample `{sample}` {i+1}/{len(files)} ({(i + 1) / len(files) * 100:.03f}%) in {elapsed:.03f}s.",
                f"Average loading time: {sum_loading_time / (i + 1):.03f}s",
                f"Average processing time: {sum_time / (i + 1):.03f}s",
                sep="\n",
                flush=True
            )

        t0 = time.time()
        col_buffer = [c for c in col_buffer if c is not None]
        if col_buffer:
            merged_tab = pd.concat(
                ([merged_tab] if merged_tab is not None else []) + col_buffer,
                axis=1
            )
            elapsed = time.time() - t0
            sum_time += elapsed
            print(f"Emptied column buffer in {elapsed:0.3f}s.", flush=True)

            # merged_tab = merged_tab.merge(column, left_index=True, right_index=True, how="outer")
            # merged_tab.rename(columns={colname: sample}, inplace=True)
            # merged_tab[sample]["unannotated"] = src_tab["uniq_raw"].get("unannotated", "NA")
            # merged_tab.loc["unannotated", sample] = src_tab["uniq_raw"].get("unannotated", "NA")
        elapsed = time.time() - t00
        print(f"Finished merging in {elapsed:0.3f}s.", flush=True)
        sum_time += elapsed

        t0 = time.time()
        merged_tab.columns = [sample for sample, _ in files]
        merged_tab.to_csv(table_file, sep="\t", na_rep="NA", index_label=feature_label)
        elapsed = time.time() - t0
        sum_time += elapsed
        print(
            f"Finished writing table in {elapsed:0.3f}s.",
            f"Average loading time: {sum_loading_time / len(files):.03f}s",
            f"Average processing time: {sum_time / len(files):.03f}s",
            sep="\n",
            flush=True,
        )

    def _collate_aln_stats(self, files):
        table_file = f"{self.prefix}.aln_stats.txt.gz"
        index = ("Total", "Passed", "Seqid", "Length")
        merged_tab = pd.DataFrame(index=index)
        for sample, fn in files:
            src_tab = pd.read_csv(fn, sep="\t", index_col=0, header=None)
            merged_tab = merged_tab.merge(src_tab, left_index=True, right_index=True, how="outer")
            merged_tab.rename(columns={1: sample}, inplace=True)
        merged_tab.to_csv(table_file, sep="\t", na_rep="NA", index_label="stat")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("count_dir", type=str)
    ap.add_argument("--out_prefix", "-o", type=str, default="./collated")
    ap.add_argument("--recursive", "-r", action="store_true")
    ap.add_argument("--column", "-c", type=str, choices=("uniq_raw", "uniq_lnorm", "uniq_scaled", "uniq_rpkm", "combined_raw", "combined_lnorm", "combined_scaled", "combined_rpkm"), default="uniq_raw")
    ap.add_argument("--index_file", type=str)
    args = ap.parse_args()

    outdir = os.path.dirname(args.out_prefix)
    if outdir and outdir != ".":
        pathlib.Path(outdir).mkdir(exist_ok=True, parents=True)

    FeatureCountCollator(args.count_dir, args.out_prefix, args.column, recursive=args.recursive).collate(index_file=args.index_file)


if __name__ == "__main__":
    main()

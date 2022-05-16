# pylint: disable=C0103, C0301
""" module docstring """

import argparse
import csv
import gzip
import os
import pathlib


def split_table(f, output_dir):
    pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)

    def write_subtable(instream, outstream):
        with outstream:
            for row in instream:
                if row[0][0] == "#":
                    return row[0][1:].strip()
                print(*row, sep="\t", file=outstream)
            return None

    with (gzip.open if f.endswith(".gz") else open)(f, "rt", encoding="utf-8") as _in:
        table_r = csv.reader(_in, delimiter="\t")
        header = next(table_r)
        unannotated = next(table_r)

        out_f = None
        fname = os.path.join(
            output_dir,
            os.path.basename(f).replace(".feature_counts", "").replace(".gz", "").replace(".txt", "")
        )

        category = None
        for row in table_r:
            write_row = False
            if row:
                if row[0][0] == "#" or category is not None:
                    if category is None:
                        category = row[0][1:].strip()
                    else:
                        write_row = True

                out_f = open(f"{fname}.{category}.txt", "w", encoding="utf-8")
                with out_f:
                    print(*header, sep="\t", file=out_f)
                    print(*unannotated, sep="\t", file=out_f)
                    if write_row:
                        print(*row, sep="\t", file=out_f)
                    category = write_subtable(table_r, out_f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table_input")
    ap.add_argument("--output_dir", "-o", default="split_output", type=str)
    args = ap.parse_args()

    split_table(args.table_input, args.output_dir)


if __name__ == "__main__":
    main()

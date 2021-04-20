import gzip
import sys
import argparse
import os

class GffIndexer:
    def __init__(self, gff, overwrite=False):
        self._refs = dict()
        offset = 0
        cur_ref = None
        index_fn = gff + ".index"
        if os.path.exists(index_fn):
            if not overwrite:
                raise FileExistsError("Index {fn} already exists. Please use -f/--force option to overwrite.".format(fn=index_fn))
            print("--force parameter is set: overwriting existing index {fn}.".format(fn=index_fn))

        with open(gff, "rt") as f, open(index_fn, "wt") as index_out:
            for line in f:
                if not line.startswith("#"):
                    ref = line.split("\t")[0]
                    if ref != cur_ref:
                        #print(self._refs)
                        if cur_ref is not None:
                            self._refs[cur_ref][-1][1] = offset - self._refs[cur_ref][-1][0]
                            print(cur_ref, *self._refs[cur_ref][-1], sep="\t", flush=True, file=index_out)
                        cur_ref = ref
                        self._refs.setdefault(cur_ref, list()).append([offset, 0])
                offset += len(line)

            self._refs[cur_ref][-1][1] = offset - self._refs[cur_ref][-1][0]
            print(cur_ref, *self._refs[cur_ref][-1], sep="\t", flush=True, file=index_out)

    def get_index(self, out=sys.stdout):
        for r, p in self._refs.items():
            for pp in p:
                print(r, *pp, file=out, sep="\t", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gff_file", type=str)
    ap.add_argument("--force", "-f", action="store_true", help="Force overwrite of existing index file (input_file.index)")
    args = ap.parse_args()

    GffIndexer(args.gff_file, overwrite=args.force)

if __name__ == "__main__":
    main()

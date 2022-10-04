""" module docstring """

import gzip

from itertools import chain

import numpy as np


class CoverageCounter(dict):
    def __init__(self):
        dict.__init__(self)

    # pylint: disable=R0913
    def update_coverage(self, rid, start, end, uniq_counts, ambig_counts, annotation):
        region_length = end - start + 1
        cov = self.setdefault(
            (rid, start, end),
            {
                "uniq_coverage": np.zeros(region_length),
                "combined_coverage": np.zeros(region_length),
                "annotation": annotation
            }
        )

        for rstart, rend, rcount in chain(*uniq_counts):
            cov["uniq_coverage"][rstart:rend] += rcount
            cov["combined_coverage"][rstart:rend] += rcount
        for rstart, rend, rcount in chain(*ambig_counts):
            cov["combined_coverage"][rstart:rend] += rcount

    def dump(self, prefix):
        with gzip.open(f"{prefix}.{self.__class__.__name__}.txt.gz", "wt") as _out:
            # pylint: disable=C0103,C0301
            for k, v in self.items():
                print(*k, ":", file=_out)
                print(
                    v["annotation"],
                    v["uniq_coverage"].mean() if v["uniq_coverage"] is not None else "NA",
                    (v["uniq_coverage"][v["uniq_coverage"] > 0]).mean() if v["uniq_coverage"] is not None else "NA",
                    v["combined_coverage"].mean() if v["combined_coverage"] is not None else "NA",
                    (v["combined_coverage"][v["uniq_coverage"] > 0]).mean() if v["combined_coverage"] is not None else "NA",
                    file=_out
                )

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
            istart, iend = rstart - start + 1, rend - start + 1  # [start, end)
            cov["uniq_coverage"][istart:iend] += rcount
            cov["combined_coverage"][istart:iend] += rcount
        for rstart, rend, rcount in chain(*ambig_counts):
            istart, iend = rstart - start + 1, rend - start + 1
            cov["combined_coverage"][istart:iend] += rcount

    def dump(self, prefix):
        with gzip.open(f"{prefix}.{self.__class__.__name__}.txt.gz", "wt") as _out:
            # pylint: disable=C0103,C0301
            for k, v in self.items():
                print(*k, ":", file=_out)
                avg_uniq_cov = v["uniq_coverage"].mean()
                avg_uniq_cov_nonzero = (v["uniq_coverage"][v["uniq_coverage"] > 0]).mean() if v["uniq_coverage"].any() else "NA"
                avg_comb_cov = v["combined_coverage"].mean()
                avg_comb_cov_nonzero = (v["combined_coverage"][v["combined_coverage"] > 0]).mean() if v["combined_coverage"].any() else "NA"
                print(
                    v["annotation"],
                    avg_uniq_cov, avg_uniq_cov_nonzero,
                    avg_comb_cov, avg_comb_cov_nonzero,
                    file=_out
                )

""" module docstring """

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

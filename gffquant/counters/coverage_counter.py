""" module docstring """

from itertools import chain

import numpy as np


class CoverageCounter(dict):
    def __init__(self):
        dict.__init__(self)

    # pylint: disable=R0913
    def update_coverage(self, rid, start, end, uniq_counts, ambig_counts, annotation):
        cov = self.setdefault(
            (rid, start, end),
            {
                "coverage": np.zeros(end - start + 1),
                "annotation": annotation
            }
        )
        for rstart, rend, rcount in chain(*uniq_counts, *ambig_counts):
            cov["coverage"][rstart:rend] += rcount

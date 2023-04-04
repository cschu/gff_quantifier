# pylint: disable=W0223
# pylint: disable=C0103
# pylint: disable=W1514

"""module docstring"""

import gzip

from collections import Counter


class AlignmentCounter(Counter):
    COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled"]

    @staticmethod
    def normalise_counts(counts, feature_len, scaling_factor):
        """Returns raw, length-normalised, and scaled feature counts."""
        normalised = counts / feature_len
        scaled = normalised * scaling_factor
        return counts, normalised, scaled

    def get_increment(self, n_aln, increment):
        # 1overN = lavern. Maya <3
        return (increment / n_aln) if self.distribution_mode == "1overN" else increment

    def __init__(self, distribution_mode="uniq_only", strand_specific=False):
        Counter.__init__(self)
        self.distribution_mode = distribution_mode
        self.strand_specific = strand_specific
        self.unannotated_reads = 0

    def dump(self, prefix, refmgr):
        with gzip.open(f"{prefix}.{self.__class__.__name__}.txt.gz", "wt") as _out:
            for k, v in self.items():
                ref, reflen = refmgr.get(k[0] if isinstance(k, tuple) else k)
                print(k, ref, reflen, v, sep="\t", file=_out)

import logging

from collections import Counter

import numpy as np

from .. import DistributionMode


logger = logging.getLogger(__name__)


class AlignmentCounter:
    COUNT_HEADER_ELEMENTS = ("raw", "lnorm", "scaled")
    INITIAL_SIZE = 1000

    @staticmethod
    def normalise_counts(counts, feature_len, scaling_factor):
        """Returns raw, length-normalised, and scaled feature counts."""
        normalised = counts / feature_len
        scaled = normalised * scaling_factor
        return counts, normalised, scaled

    def get_increment(self, n_aln, increment):
        # 1overN = lavern. Maya <3
        return (increment / n_aln) if self.distribution_mode == DistributionMode.ONE_OVER_N else increment

    def __init__(self, distribution_mode=DistributionMode.ONE_OVER_N, strand_specific=False):
        self.distribution_mode = distribution_mode
        self.strand_specific = strand_specific
        self.unannotated_reads = 0

        self.index = {}
        self.counts = np.zeros(
            (AlignmentCounter.INITIAL_SIZE, 2,),
            dtype='float64',
        )
    def dump(self, prefix, refmgr):
        import gzip
        with gzip.open(f"{prefix}.{self.__class__.__name__}.txt.gz", "wt") as _out:
            for key in self:
                ref, reflen = refmgr.get(key[0] if isinstance(key, tuple) else key)
                print(key, ref, reflen, self.counts[key], sep="\t", file=_out)
            # for k, v in self.items():
            # ref, reflen = refmgr.get(k[0] if isinstance(k, tuple) else k)
            # print(k, ref, reflen, v, sep="\t", file=_out)
        ...
    def get(self, key, default_val):
        key_index = self.index.get(key)
        if key_index is None:
            return Counter()
        return Counter({key: self.counts[key_index]})
    
    def setdefault(self, key, default_val):
        ...

    def has_ambig_counts(self):
        return bool(self.counts[:, 1].sum() != 0)
    
    def __iter__(self):
        yield from self.index.keys()
    def __getitem__(self, key):
        key_index = self.index.get(key)
        if key_index is None:
            return 0.0
        return self.counts[key_index]
    def __setitem__(self, key, value):
        key_index = self.index.get(key)
        if key_index is not None:
            self.counts[key_index] = value
        else:
            raise KeyError(f"{key=} not found.")
        
    def update_counts(self, count_stream, increment=1, ambiguous_counts=False):
        contributed_counts = 0
        for hits, aln_count in count_stream:
            hit = hits[0]
            inc = (
                (
                    self.get_increment(aln_count, increment),
                    increment,
                )
            )[aln_count == 1]
            key = (
                (
                    (hit.rid, hit.rev_strand),
                    hit.rid
                )
            )[self.strand_specific]

            key_index = self.index.get(key)
            if key_index is None:
                nrows = self.counts.shape[0]
                if len(self.index) == nrows:
                    self.counts = np.pad(
                        self.counts,
                        ((0, AlignmentCounter.INITIAL_SIZE), (0, 0),),
                    )
                # key_index = self.index.setdefault(key, len(self.index))
                key_index = self.index[key] = len(self.index)
            self.counts[key_index][int(ambiguous_counts)] += inc
            contributed_counts += inc

        return contributed_counts
    
    def transform(self, refmgr):
        # transform 2-column uniq/ambig count matrix
        # into 4 columns
        # uniq_raw, combined_raw, uniq_lnorm, combined_lnorm

        # obtain gene lengths
        lengths = np.array(
            tuple(
                (refmgr.get(key[0] if isinstance(key, tuple) else key))[1]
                for key in self.index
            )
        )
        logger.info("LENGTHS ARRAY = %s", lengths.shape)
        logger.info("INDEX SIZE = %s", len(self.index))

        # remove the un-indexed rows
        self.counts = self.counts[0:len(self.index), :]

        # calculate combined_raw
        self.counts[:, 1:2] += self.counts[:, 0:1]

        # duplicate the raw counts
        self.counts = np.concatenate(
            (self.counts, self.counts,),
            axis=1,
        )

        # length-normalise the lnorm columns
        self.counts[:, 2:4] /= lengths[:, None]

        # return count sums
        return self.counts.sum(axis=0)



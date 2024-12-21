import logging

from collections import Counter

import numpy as np

from .. import DistributionMode


logger = logging.getLogger(__name__)


class AlignmentCounter:
    COUNT_HEADER_ELEMENTS = ("raw", "lnorm", "scaled")
    INITIAL_SIZE = 1000
    PLUS_STRAND, MINUS_STRAND = False, True

    @staticmethod
    def normalise_counts(counts, feature_len, scaling_factor):
        """Returns raw, length-normalised, and scaled feature counts."""
        normalised = counts / feature_len
        scaled = normalised * scaling_factor
        return counts, normalised, scaled

    def get_increment(self, n_aln, increment):
        # 1overN = lavern. Maya <3
        return (increment / n_aln) if self.distribution_mode == DistributionMode.ONE_OVER_N else increment
    
    def toggle_single_read_handling(self, unmarked_orphans):
        # precalculate count-increment for single-end, paired-end reads
        # for mixed input (i.e., paired-end data with single-end reads = orphans from preprocessing),
        # properly attribute fractional counts to the orphans
        # Increments:
        # alignment from single end library read: 1
        # alignment from paired-end library read: 0.5 / mate (pe_count = 1) or 1 / mate (pe_count = 2)
        # alignment from paired-end library orphan: 0.5 (pe_count = 1) or 1 (pe_count = 2)

        # old code:
        # increment = 1 if (not pair or self.paired_end_count == 2) else 0.5

        # if pair:
        #     increment = 1 if self.paired_end_count == 2 else 0.5
        # else:
        #     increment = 0.5 if self.unmarked_orphans else 1
        self.increments = (
            (self.paired_end_count / 2.0) if unmarked_orphans else 1.0,
            self.paired_end_count / 2.0,
        )

    def __init__(
        self,
        distribution_mode=DistributionMode.ONE_OVER_N,
        strand_specific=False,
        paired_end_count=1,
    ):
        self.distribution_mode = distribution_mode
        self.strand_specific = strand_specific
        self.paired_end_count = paired_end_count
        self.increments = (1.0, 1.0,)
        self.increments_auto_detect = (1.0, self.paired_end_count / 2.0,)
        self.unannotated_reads = 0

        self.index = {}
        self.counts = np.zeros(
            (AlignmentCounter.INITIAL_SIZE, 2,),
            dtype='float64',
        )
    def dump(self, prefix, refmgr):
        import gzip
        with gzip.open(f"{prefix}.{self.__class__.__name__}.txt.gz", "wt") as _out:
            for key, key_index in self.index.items():
                ref, reflen = refmgr.get(key[0] if isinstance(key, tuple) else key)
                print(key, ref, reflen, self.counts[key_index], sep="\t", file=_out)
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
    
    def update(self, count_stream, ambiguous_counts=False, pair=False, pe_library=None,):
        if pe_library is not None:
            # this is the case when the alignment has a read group tag
            # if pe_library is True (RG tag '2') -> take paired-end increment (also for orphans)
            # else (RG tag '1') -> take single-end increment
            increment = self.increments_auto_detect[pe_library]
        else:
            # if the alignment has no (appropriate) read group tag
            # use the paired-end information instead
            # if orphan reads are present in the input sam/bam,
            # the flag `--unmarked_orphans` should be set
            # otherwise orphan reads will be assigned a count of 1.
            increment = self.increments[pair]

        contributed_counts = self.update_counts(count_stream, increment=increment, ambiguous_counts=ambiguous_counts,)

        return contributed_counts
    
    def get_unannotated_reads(self):
        return self.unannotated_reads
    
    def get_counts(self, seqid, strand_specific=False):
        if strand_specific:
            raise NotImplementedError()
        counts = self[seqid]
        return np.array((counts[0], counts[2], counts[1], counts[3]))
    
    def get_all_regions(self):
        yield from self 
        
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
                    hit.rid,
                    (hit.rid, hit.rev_strand),
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



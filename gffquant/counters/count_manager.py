"""module docstring"""

from collections import Counter

from .region_counter import UniqueRegionCounter, AmbiguousRegionCounter
from .seq_counter import UniqueSeqCounter, AmbiguousSeqCounter


# pylint: disable=R0902
class CountManager:
    # this may be counter-intuitive
    # but originates from the samflags 0x10, 0x20,
    # which also identify the reverse-strandness of the read
    # and not the forward-strandness
    PLUS_STRAND, MINUS_STRAND = False, True

    def __init__(
        # pylint: disable=W0613,R0913
        self,
        distribution_mode="1overN",
        region_counts=True,
        strand_specific=False,
        calc_coverage=False,
        paired_end_count=1,
        unmarked_orphans=False,
    ):
        self.distribution_mode = distribution_mode
        self.strand_specific = strand_specific
        self.calc_coverage = calc_coverage
        # precalculate count-increment for single-end, paired-end reads
        # for mixed input (i.e., paired-end data with single-end reads = orphans from preprocessing),
        # properly attribute fractional counts to the orphans
        # Increments:
        # alignment from single end library read: 1
        # alignment from paired-end library read: 0.5 / mate (pe_count = 1) or 1 / mate (pe_count = 2)
        # alignment from paired-end library orphan: 0.5 (pe_count = 1) or 1 (pe_count = 2)
        self.increments = [
            (paired_end_count / 2.0) if unmarked_orphans else 1.0,
            paired_end_count / 2.0
        ]

        self.uniq_seqcounts, self.ambig_seqcounts = None, None
        self.uniq_regioncounts, self.ambig_regioncounts = None, None

        if region_counts:
            self.uniq_regioncounts = UniqueRegionCounter(strand_specific=strand_specific, calc_coverage=calc_coverage)
            self.ambig_regioncounts = AmbiguousRegionCounter(
                strand_specific=strand_specific,
                distribution_mode=distribution_mode,
                calc_coverage=calc_coverage
            )

        else:
            self.uniq_seqcounts = UniqueSeqCounter(strand_specific=strand_specific)
            self.ambig_seqcounts = AmbiguousSeqCounter(
                strand_specific=strand_specific,
                distribution_mode=distribution_mode
            )

    def has_ambig_counts(self):
        return self.ambig_regioncounts or self.ambig_seqcounts

    def update_counts(self, count_stream, ambiguous_counts=False, pair=False):
        seq_counter, region_counter = (
            (self.uniq_seqcounts, self.uniq_regioncounts)
            if not ambiguous_counts
            else (self.ambig_seqcounts, self.ambig_regioncounts)
        )

        increment = self.increments[pair]

        # increment = 1 if (not pair or self.paired_end_count == 2) else 0.5

        # if pair:
        #     increment = 1 if self.paired_end_count == 2 else 0.5
        # else:
        #     increment = 0.5 if self.unmarked_orphans else 1

        if seq_counter is not None:
            seq_counter.update_counts(count_stream, increment=increment)
        elif region_counter is not None:
            region_counter.update_counts(count_stream, increment=increment)

    def dump_raw_counters(self, prefix, bam):
        if self.uniq_seqcounts is not None:
            self.uniq_seqcounts.dump(prefix, bam)
        if self.ambig_seqcounts is not None:
            self.ambig_seqcounts.dump(prefix, bam)
        if self.uniq_regioncounts is not None:
            self.uniq_regioncounts.dump(prefix, bam)
        if self.ambig_regioncounts is not None:
            self.ambig_regioncounts.dump(prefix, bam)

    def get_unannotated_reads(self):
        unannotated_reads = 0

        if self.uniq_regioncounts is not None:
            unannotated_reads += self.uniq_regioncounts.unannotated_reads
        if self.ambig_regioncounts is not None:
            unannotated_reads += self.ambig_regioncounts.unannotated_reads
        if self.uniq_seqcounts is not None:
            unannotated_reads += self.uniq_seqcounts.unannotated_reads
        if self.ambig_seqcounts is not None:
            unannotated_reads += self.ambig_seqcounts.unannotated_reads

        return unannotated_reads

    def get_counts(self, seqid, region_counts=False, strand_specific=False):
        if region_counts:
            rid, seqid = seqid[0], seqid[1:]
            uniq_counter = self.uniq_regioncounts.get(rid, {} if self.calc_coverage else Counter())
            ambig_counter = self.ambig_regioncounts.get(rid, {} if self.calc_coverage else Counter())

            # pylint: disable=R1720
            if strand_specific:
                raise NotImplementedError
            elif self.calc_coverage:
                return [uniq_counter.get(seqid, [])], [ambig_counter.get(seqid, [])]
            else:
                return [uniq_counter[seqid]], [ambig_counter[seqid]]

        else:
            uniq_counter, ambig_counter = self.uniq_seqcounts, self.ambig_seqcounts

            if strand_specific:
                uniq_counts, ambig_counts = [0.0, 0.0], [0.0, 0.0]
                uniq_counts[seqid[1]] = uniq_counter[seqid]
                ambig_counts[seqid[1]] = ambig_counter[seqid]

                # rid = seqid[0] if isinstance(seqid, tuple) else seqid
                # uniq_counts = [
                #     uniq_counter[(rid, CountManager.PLUS_STRAND)],
                #     uniq_counter[(rid, CountManager.MINUS_STRAND)],
                # ]
                # ambig_counts = [
                #     ambig_counter[(rid, CountManager.PLUS_STRAND)],
                #     ambig_counter[(rid, CountManager.MINUS_STRAND)],
                # ]
            else:
                uniq_counts, ambig_counts = [uniq_counter[seqid]], [ambig_counter[seqid]]

            return uniq_counts, ambig_counts

    def get_regions(self, rid):
        return set(self.uniq_regioncounts.get(rid, set())).union(
            self.ambig_regioncounts.get(rid, set())
        )

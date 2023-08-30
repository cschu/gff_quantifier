# pylint: disable=W0223

"""module docstring"""

from collections import Counter

from .alignment_counter import AlignmentCounter


class RegionCounter(AlignmentCounter):
    """This counter class can be used in overlap mode, i.e.
    when reads are aligned against long references (e.g. contigs)
    with multiple regions of interest (features).
    """

    def __init__(self, distribution_mode="uniq_only", strand_specific=False):
        AlignmentCounter.__init__(
            self, distribution_mode=distribution_mode, strand_specific=strand_specific
        )

    # pylint: disable=R0913,W0613
    def _update_region(self, region_id, ostart, oend, rev_strand, cstart, cend, increment=1):
        overlap_id = ((ostart, oend), rev_strand) if self.strand_specific else (ostart, oend)
        self.setdefault(region_id, Counter())[overlap_id] += increment


class UniqueRegionCounter(RegionCounter):
    """This counter class can be used in overlap mode, i.e.
    when reads are aligned against long references (e.g. contigs)
    with multiple regions of interest (features).
    """

    def __init__(self, distribution_mode="uniq_only", strand_specific=False):
        RegionCounter.__init__(
            self, distribution_mode=distribution_mode, strand_specific=strand_specific,
        )

    # pylint: disable=W0613
    def update_counts(self, count_stream, increment=1, pair=False):
        """Update counter with alignments against the same reference.

        input: count_stream
        - counts: set of overlaps with the reference
        - aln_count: 1 if overlaps else 0
        - unaligned: 1 - aln_count
        (redundant input due to streamlining uniq/ambig dataflows)
        """
        for counts, aln_count, unaligned in count_stream:
            if aln_count:
                for rid, hits in counts.items():
                    for hit in hits:
                        self._update_region(
                            rid, *hit, increment=increment
                        )
            else:
                self.unannotated_reads += unaligned


class AmbiguousRegionCounter(RegionCounter):
    """This counter class can be used in overlap mode, i.e.
    when reads are aligned against long references (e.g. contigs)
    with multiple regions of interest (features).
    """

    def __init__(self, distribution_mode="uniq_only", strand_specific=False):
        RegionCounter.__init__(
            self, distribution_mode=distribution_mode, strand_specific=strand_specific,
        )

    # pylint: disable=W0613
    def update_counts(self, count_stream, pair=False, increment=1):
        """Update counter with alignments against the same reference.

        input: count_stream
        - counts: set of overlaps with the reference
        - aln_count: 1 if overlaps else 0
        - unaligned: 1 - aln_count
        (redundant input due to streamlining uniq/ambig dataflows)
        """
        for counts, aln_count, unaligned in count_stream:
            if aln_count:
                inc = self.get_increment(aln_count, increment)
                for rid, hits in counts.items():
                    for hit in hits:
                        self._update_region(
                            rid, *hit, increment=inc
                        )
            else:
                self.unannotated_reads += unaligned

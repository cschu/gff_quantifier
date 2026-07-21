# pylint: disable=W0223,R0917

"""module docstring"""

from collections import Counter

from .. import DistributionMode
from .alignment_counter import AlignmentCounter


# from count_manager.get_counts()
# if region_counts:
#     raise NotImplementedError()
#     rid, seqid = seqid[0], seqid[1:]

#     uniq_counter = self.uniq_regioncounts.get(rid, Counter())
#     ambig_counter = self.ambig_regioncounts.get(rid, Counter())

#     # pylint: disable=R1720
#     if strand_specific:
#         raise NotImplementedError
#     else:
#         return [uniq_counter[seqid]], [ambig_counter[seqid]]


class RegionCounter(AlignmentCounter):
    """This counter class can be used in overlap mode, i.e.
    when reads are aligned against long references (e.g. contigs)
    with multiple regions of interest (features).
    """

    def __init__(self, distribution_mode=DistributionMode.ONE_OVER_N, strand_specific=False):
        AlignmentCounter.__init__(
            self, distribution_mode=distribution_mode, strand_specific=strand_specific
        )

    # pylint: disable=R0913,W0613
    def _update_region(self, region_id, ostart, oend, rev_strand, cstart=None, cend=None, increment=1):
        overlap_id = ((ostart, oend), rev_strand) if self.strand_specific else (ostart, oend)
        self.setdefault(region_id, Counter())[overlap_id] += increment

    def update_counts(self, count_stream, increment=1):
        contributed_counts = 0
        for hits, aln_count in count_stream:
            inc = increment if aln_count == 1 else AlignmentCounter.get_increment(aln_count, increment, self.distribution_mode)
            for hit in hits:
                self._update_region(
                    hit.rid, hit.start, hit.end, hit.rev_strand, increment=inc,
                )
                contributed_counts += inc
        return contributed_counts

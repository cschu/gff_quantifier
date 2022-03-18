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

    def _update_region(self, region_id, overlap_id, increment=1):
        self.setdefault(region_id, Counter())[
            (overlap_id[:2], overlap_id[2])
        ] += increment

    def update_counts(self, count_stream):
        """Update counter with alignments against the same reference.

        input: count_stream
        - counts: set of overlaps with the reference
        - aln_count: 1 if overlaps else 0
        - unaligned: 1 - aln_count
        (redundant input due to streamlining uniq/ambig dataflows)
        """
        for counts, aln_count, unaligned in count_stream:
            increment = (
                (1 / aln_count) if self.distribution_mode == "1overN" else 1
            )  # 1overN = lavern. Maya <3
            for rid, hits in counts.items():
                for ostart, oend, rev_strand, _, _ in hits:
                    self._update_region(
                        rid, (ostart, oend, rev_strand), increment=increment
                    )

            self.unannotated_reads += unaligned
            yield counts, aln_count, unaligned

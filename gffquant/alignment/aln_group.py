""" module docstring """
from itertools import chain


class AlignmentGroup:
    def __init__(self, aln=None):
        self.qname = None
        self.primaries = [None, None]
        self.secondaries = [[], []]

        if aln is not None:
            self.add_alignment(aln)

    def add_alignment(self, aln):
        if self.qname is None:
            self.qname = aln.qname
        if aln.is_unique():
            self.primaries[aln.is_second()] = aln
        else:
            self.secondaries[aln.is_second()].append(aln)

    def n_align(self):
        return sum(
            (self.primaries[0] is not None),
            (self.primaries[1] is not None),
            len(self.secondaries[0]),
            len(self.secondaries[1])
        )

    def get_alignments(self):
        yield self.primaries[0]
        yield self.primaries[1]
        for aln in chain(*self.secondaries):
            yield aln

    def get_ambig_align_counts(self):
        return (
            bool(self.primaries[False]) + len(self.secondaries[False]),
            bool(self.primaries[True]) + len(self.secondaries[True])
        ) if len(self.secondaries[False]) + len(self.secondaries[True]) else (0, 0)

    def is_aligned_pair(self):
        aln1, aln2 = self.primaries
        return all(
            (
                aln1 is not None and aln2 is not None,
                aln1.rid == aln2.rid,
                aln1.is_paired() and aln2.is_paired()
            )
        )

""" module docstring """
from itertools import chain


class AlignmentGroup:
    def __init__(self, aln=None):
        self.qname = None
        self.primaries = [None, None]
        self.secondaries = [[], []]
        self.pe_library = None
        # hit counts are mode-dependent:
        # ambiguous alignments in true ambig modes require the number of alternative locations
        # (we handle first/second PE-alignments separately, hence a two-element list)
        # otherwise the hit count is always 1
        self.ambig_hit_counts = [0, 0]

        if aln is not None:
            self.add_alignment(aln)

    def add_alignment(self, aln):
        if self.qname is None:
            self.qname = aln.qname
        if aln.is_primary():
            self.primaries[aln.is_second()] = aln
        else:
            self.secondaries[aln.is_second()].append(aln)

        if aln.hits:
            self.ambig_hit_counts[aln.is_second()] += 1

        # alignments generated by a native alignment runner or from a nevermore upstream process
        # will have paired-end information in the RG-tag

        # if the incoming alignment has no RG-tag but the aln-group has pe_library set,
        # which can only happen via a previously added read with read group
        # -> conflict
        if aln.read_group is None and self.pe_library is not None:
            raise ValueError(f"Alignment {str(aln)} has no read group information. Expected: {self.pe_library}")

        # now check, if the alignment has a valid RG-tag
        # then compare it with the pe_library of the aln_group
        # set pe_library if not set else check for conflict
        if aln.read_group is not None and aln.read_group in (1, 2):
            pe_library = aln.read_group == 2
            if self.pe_library is not None and pe_library != self.pe_library:
                raise ValueError(f"Conflicting read group information found in {str(aln)}.")
            self.pe_library = pe_library

    def n_align(self):
        return sum((
            (self.primaries[0] is not None),
            (self.primaries[1] is not None),
            len(self.secondaries[0]),
            len(self.secondaries[1])
        ))

    def get_alignments(self):
        for aln in chain(self.primaries, *self.secondaries):
            if aln is not None:
                yield aln

    def get_all_hits(self, as_ambiguous=False):
        is_ambiguous = as_ambiguous and self.is_ambiguous()
        for aln in self.get_alignments():
            if aln.hits:
                n_aln = (1, self.ambig_hit_counts[aln.is_second()])[is_ambiguous]
                yield aln.hits, n_aln

    def get_ambig_align_counts(self):
        return (
            bool(self.primaries[False]) + len(self.secondaries[False]),
            bool(self.primaries[True]) + len(self.secondaries[True])
        ) if len(self.secondaries[False]) + len(self.secondaries[True]) else (0, 0)

    def is_ambiguous(self):
        return any(self.secondaries)

    def is_aligned_pair(self):
        aln1, aln2 = self.primaries

        try:
            return aln1.rid == aln2.rid and aln1.is_paired() and aln2.is_paired()
        except AttributeError:
            return False

    def is_paired(self):
        return any(aln.is_paired() for aln in self.get_alignments())

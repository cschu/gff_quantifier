""" module docstring """

from ..bamreader import BamFile, SamFlags


class PairedEndAlignmentCache(dict):
    def __init__(self, ambig_alignments=False):
        dict.__init__(self)
        self.ambig_alignments = ambig_alignments

    def empty_cache(self):
        for qname, alignment in self.items():
            n_aln = len(alignment)

            if n_aln == 1:
                #  unpaired
                start, end, flag = alignment[0][1:]
                rev_strand = SamFlags.is_reverse_strand(flag)
            elif n_aln == 2:
                #  paired
                start, end = BamFile.calculate_fragment_borders(
                    *alignment[0][1:-1], *alignment[1][1:-1]
                )
                # pylint: disable=W0511
                rev_strand = (
                    None  #  TODO: add strand-specific handling by RNAseq protocol
                )
            else:
                raise ValueError(
                    f"{n_aln} primary alignments detected for read {qname}."
                )

            yield alignment[0][0], start, end, rev_strand

        self.clear()

    def process_alignment(self, alignment, merge_pair=False):
        """ returns the mate from the cache or merges the pair """
        # need to keep track of read pairs to avoid dual counts

        start, end, rev_strand = None, None, None

        # check if the mate has already been seen
        mates = self.setdefault(alignment.qname, [])
        if mates:
            if len(mates) > 1:
                str_mates = '\n'.join(mates)
                raise ValueError(
                    f"Alignment {alignment.qname} has more than two mates.\n{str_mates}"
                )

            rid, start, end, flag = mates[0]

            if alignment.rnext != rid:
                raise ValueError(
                    f"Alignment {alignment.qname} seems to be corrupted: "
                    f"{str(alignment)} {rid}."
                )

            # if requested, merge the pair
            if merge_pair:
                start, end = BamFile.calculate_fragment_borders(
                    alignment.start, alignment.end, *mates[0][1:-1]
                )
            else:
                rev_strand = SamFlags.is_reverse_strand(flag)

            # and remove the pair from the cache
            del self[alignment.qname]
        else:
            # otherwise cache the first encountered mate and advance to the next read
            mates.append(
                (alignment.rid, alignment.start, alignment.end, alignment.flag)
            )

        return start, end, rev_strand

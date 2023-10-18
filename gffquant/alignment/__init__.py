# flake8: noqa

""" module docstring """

from dataclasses import dataclass

from .aln_group import AlignmentGroup
from .pysam_alignment_processor import AlignmentProcessor
from .samflags import SamFlags
from .cigarops import CigarOps


@dataclass
class EncodedAlignment:
    """ class to reduce alignment data """
    qname: int = None
    flag: int = None
    rid: int = None
    start: int = None
    end: int = None
    read_group: int = None

    @classmethod
    def from_pysam_alignment(cls, pysam_aln, read_id, rid):
        try:
            rg_tag = pysam_aln.get_tag("RG")
        except KeyError:
            rg_tag = None

        return cls(
            read_id,
            pysam_aln.flag,
            rid,
            pysam_aln.pos,
            CigarOps.calculate_coordinates(
                pysam_aln.pos,
                [(y, x) for x, y in pysam_aln.cigar]
            ),
            rg_tag,
        )

    # various flag checks -- only needed ones are implemented
    def is_primary(self):
        """ is this flagged as primary alignment? """
        return not bool(self.flag & SamFlags.SECONDARY_ALIGNMENT)

    def is_second(self):
        """ is this flagged as second in pair? """
        return bool(self.flag & SamFlags.SECOND_IN_PAIR == SamFlags.SECOND_IN_PAIR)

    def is_paired(self):
        """ is this flagged as paired? """
        return bool(self.flag & SamFlags.PAIRED)

    def is_reverse(self):
        """ is this flagged as reversed? """
        return bool(self.flag & SamFlags.REVERSE)

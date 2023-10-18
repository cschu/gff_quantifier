# flake8: noqa

""" module docstring """

from dataclasses import dataclass

from .aln_group import AlignmentGroup
from .pysam_alignment_processor import AlignmentProcessor
from .samflags import SamFlags
from .cigarops import CigarOps


@dataclass
class Alignment:
    read_id: int = None
    flag: int = None
    reference_id: int = None
    start: int = None
    end: int = None
    mapq: int = None
    length: int = None
    counts: int = None

    @classmethod
    def from_pysam_alignment(cls, pysam_aln, read_id, rid):
        return cls(
            read_id,
            pysam_aln.flag,
            rid,
            pysam_aln.pos,
            CigarOps.calculate_coordinates(
                pysam_aln.pos,
                [(y, x) for x, y in pysam_aln.cigar]
            ),
            pysam_aln.mapq,
            pysam_aln.alen,
            dict(pysam_aln.tags).get("RG")
        )

# pylint: disable=R0913

""" module docstring """

from enum import Enum, auto, unique


__version__ = "2.16.2"
__tool__ = "gffquant"


@unique
class RunMode(Enum):
    GENE = (auto(), False, False)
    DOMAIN = (auto(), True, False)
    SMALL_GENOME = (auto(), True, True)

    def __init__(self, num, overlap_required, report_unannotated):
        self.num = num
        self.overlap_required = overlap_required
        self.report_unannotated = report_unannotated

    @classmethod
    def parse(cls, string):
        string = string.upper().replace(" ", "_")
        return cls.__members__.get(
            string,
            cls.__members__.get(f"{string[:-1]}")
        )


@unique
class DistributionMode(Enum):
    ONE_OVER_N = (auto(), "1overN", True, True, True)
    ALL_ONE = (auto(), "all1", True, True, False)
    UNIQUE_ONLY = (auto(), "unique_only", False, False, False)
    PRIMARY_ONLY = (auto(), "primary_only", True, False, False)

    def __init__(self, num, alias, allow_ambiguous, allow_secondary, require_ambig_tracking):
        self.num = num
        self.alias = alias
        self.allow_ambiguous = allow_ambiguous
        self.allow_secondary = allow_secondary
        # ambig tracking: all alignments of a read need to be processed together
        self.require_ambig_tracking = require_ambig_tracking

    @classmethod
    def parse(cls, string):
        for member in cls.__members__.values():
            if member.alias == string:
                return member
        return None

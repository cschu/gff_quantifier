""" module docstring """

from enum import Enum, auto, unique


__version__ = "2.15.0"
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

# pylint: disable=C0103,R0902,R0913,W2301

""" module docstring """

from dataclasses import dataclass, field


@dataclass
class DefaultDatabaseInputFormat:
    """ Default database input format. """
    offsets: tuple = field(default=(0, 0))
    columns: tuple = field(default=(0, 1, 2, 3))
    separator: str = "\t"


@dataclass
class BedDatabaseInputFormat(DefaultDatabaseInputFormat):
    """ BED database input format. """
    # we store everything as 1-based, closed intervals internally
    # bed coords coming in as [x,y)_0 -> [x+1, y]_1
    offsets: tuple = field(default=(1, 0))


@dataclass
class HmmerDatabaseInputFormat(DefaultDatabaseInputFormat):
    """ HMMer database input format. """
    columns: tuple = field(default=(0, 1, 2, 4))
    separator: str = ","


DB_SETTINGS_SELECTION = {
    "default": DefaultDatabaseInputFormat,
    "bed": BedDatabaseInputFormat,
    "hmmer": HmmerDatabaseInputFormat,
}

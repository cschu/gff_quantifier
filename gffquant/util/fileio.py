""" module docstring """

import gzip


def get_open_function(f):
    """ Returns a file open function corresponding to gzip-compression status. """
    gz_magic = b"\x1f\x8b\x08"
    # pylint: disable=R1732,W0511
    gzipped = open(f, "rb").read(3).startswith(gz_magic)
    return gzip.open if gzipped else open

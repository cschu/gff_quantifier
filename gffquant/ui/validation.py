# pylint: disable=C0103

""" validation """
import logging
import os


logger = logging.getLogger(__name__)


def check_minimap2_index(f):
    """ docstring """
    return os.path.isfile(f) and f.endswith(".mmi")


def check_bwa_index(prefix):
    """ docstring """
    suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")
    return all(os.path.isfile(prefix + suffix) for suffix in suffixes)


def check_input_reads(fwd_reads=None, rev_reads=None, single_reads=None, orphan_reads=None):
    """ docstring """

    all_readsets = []

    if fwd_reads and rev_reads:
        if len(fwd_reads) == len(rev_reads):
            all_readsets += zip(
                (["paired"] * len(fwd_reads)),
                fwd_reads, rev_reads
            )
        else:
            raise ValueError(
                f"Found different numbers of forward/R1 {len(fwd_reads)} "
                f"and reverse/R2 {len(rev_reads)} reads."
            )
    elif fwd_reads:
        logger.warning(
            "Found -1 forward/R1 reads but no -2 reverse/R2 reads. "
            "Treating these as single-end reads."
        )
        all_readsets += zip((["single"] * len(fwd_reads)), fwd_reads)
    elif rev_reads:
        logger.warning(
            "Found -2 reverse/R2 reads but no -1 forward/R1 reads. "
            "Treating these as single-end reads."
        )
        all_readsets += zip((["single"] * len(rev_reads)), rev_reads)

    if single_reads:
        all_readsets += zip((["single"] * len(single_reads)), single_reads)
    if orphan_reads:
        all_readsets += zip((["orphan"] * len(orphan_reads)), orphan_reads)

    if not all_readsets:
        raise ValueError("No input reads specified.")

    for _, *reads in all_readsets:
        for r in reads:
            if not os.path.isfile(r):
                raise ValueError(f"{r} does not seem to be a valid read file.")

    return all_readsets

# pylint: disable=C0103,C0301

""" module docstring """

import logging
import os
import pathlib
import subprocess
import sys

# pylint: disable=W0611
from .db.db_import import DomainBedDatabaseImporter
from .handle_args import handle_args
from .profilers import GeneQuantifier, RegionQuantifier
from .runners.alignment_runner import BwaMemRunner, Minimap2Runner
from . import __version__


logger = logging.getLogger(__name__)


def check_input_reads(fwd_reads=None, rev_reads=None, single_reads=None, orphan_reads=None):
    """ docstring """
    # fwd_reads = fwd if fwd else None
    # rev_reads = rev if rev else None
    # single_reads = singles if singles else None
    # orphan_reads = orphans if orphans else None

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


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
    logger.info("Command: %s %s", os.path.basename(sys.argv[0]), " ".join(sys.argv[1:]))

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
        )

    kwargs = {}
    annotation_db = args.annotation_db
    if args.mode in ("gene", "genes"):
        Quantifier = GeneQuantifier
    else:
        Quantifier, kwargs["reference_type"] = RegionQuantifier, args.mode
        if args.mode == "domain":
            annotation_db = DomainBedDatabaseImporter(
                logger, args.annotation_db, single_category="feature"
            )
            logger.info("Finished loading database.")

    input_data = check_input_reads(
        fwd_reads=args.reads1, rev_reads=args.reads2,
        single_reads=args.singles, orphan_reads=args.orphans,
    )

    profiler = Quantifier(
        db=annotation_db,
        out_prefix=args.out_prefix,
        ambig_mode=args.ambig_mode,
        strand_specific=args.strand_specific,
        paired_end_count=args.paired_end_count,
        **kwargs,
    )

    if args.input_type == "fastq":

        Aln_runner = {
            "bwa": BwaMemRunner,
            "minimap2": Minimap2Runner,
        }.get(args.aligner)

        if Aln_runner is None:
            raise ValueError(f"Aligner `{args.aligner}` is not supported.")

        aln_runner = Aln_runner(
            args.cpus_for_alignment,
            args.reference,
            sample_id=os.path.basename(args.out_prefix),
        )        

        for input_type, *reads in input_data:

            logger.info("Running %s alignment: %s", input_type, ",".join(reads))
            stream = aln_runner.run(reads, logger, single_end_reads=input_type != "orphan", min_identity=args.min_identity, min_seqlen=args.min_seqlen, alignment_file=args.keep_alignment_file)

            profiler.count_alignments(
                stream, aln_format="sam", min_identity=args.min_identity, min_seqlen=args.min_seqlen,
            )

    else :
        profiler.count_alignments(
            args.input_files,
            aln_format=args.input_type,
            min_identity=args.min_identity,
            min_seqlen=args.min_seqlen,
            external_readcounts=args.import_readcounts,
            unmarked_orphans=args.unmarked_orphans,
        )

    profiler.finalise(
        report_category=True,
        report_unannotated=args.mode == "genome",
        dump_counters=args.debug,
        restrict_reports=args.restrict_metrics,
    )


if __name__ == "__main__":
    main()

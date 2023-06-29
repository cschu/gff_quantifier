# pylint: disable=C0103,C0301

""" module docstring """

import logging
import os
import pathlib
import subprocess
import sys

# pylint: disable=W0611
from gffquant.db.db_import import DomainBedDatabaseImporter
from gffquant.profilers import GeneQuantifier, RegionQuantifier
from . import __version__
from .handle_args import handle_args


logger = logging.getLogger(__name__)


# pylint: disable=R0913
def run_alignment(
    profiler,
    input_files,
    bwa_index,
    cpus_for_alignment=1,
    min_identity=None,
    min_seqlen=None,
    single_end_reads=False,
    blocksize=10000000,
    sample_id="sample_x",
):
    """ docstring """

    # def read_group_id = (sample.library == "paired") ? ((sample.is_paired) ? 2 : 2) : 1
    # def read_group = "'@RG\\tID:${read_group_id}\\tSM:${sample.id}'"
    # -R ${read_group}

    read_group = f"'@RG\tID:{1 if single_end_reads else 2}\tSM:{sample_id}'"

    commands = (
        f"cat {input_files} ",
        f"bwa mem -p -v 1 -a -t {cpus_for_alignment} -R {read_group}"
        f"-K {blocksize} {bwa_index} -",
    )
    

    commands = " | ".join(commands)

    logger.info("Used command: %s", commands)

    try:
        with subprocess.Popen(
            commands, shell=True, stdout=subprocess.PIPE
        ) as read_processing_proc:
            profiler.count_alignments(
                read_processing_proc.stdout,
                aln_format="sam",
                min_identity=min_identity,
                min_seqlen=min_seqlen,                
            )
    except Exception as err:
        logger.error("Caught some exception:")
        logger.error("%s", err)
        raise Exception from err


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
        quantifier = GeneQuantifier
    else:
        quantifier, kwargs["reference_type"] = RegionQuantifier, args.mode
        if args.mode == "domain":
            annotation_db = DomainBedDatabaseImporter(
                logger, args.annotation_db, single_category="cazy"
            )
            logger.info("Finished loading database.")

    profiler = quantifier(
        db=annotation_db,
        out_prefix=args.out_prefix,
        ambig_mode=args.ambig_mode,
        strand_specific=args.strand_specific,
        paired_end_count=args.paired_end_count,
        **kwargs,
    )

    if args.input_type == "fastq":
        run_alignment(
            profiler,
            args.input_files,
            args.reference,
            args.cpus_for_alignment,
            args.min_identity,
            args.min_seqlen,
            single_end_reads=not args.unmarked_orphans,
            sample_id=os.path.basename(args.output_prefix),
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

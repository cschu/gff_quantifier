# pylint: disable=C0103,C0301

""" module docstring """

import logging
import os
import sys

# pylint: disable=W0611
from .db.importers import SmallDatabaseImporter, SmallGenomeDatabaseImporter
from .handle_args import handle_args
from .profilers import GeneQuantifier, RegionQuantifier, FeatureQuantifier
from .runners.alignment_runner import BwaMemRunner, Minimap2Runner
from . import __version__, RunMode


logger = logging.getLogger(__name__)


def stream_alignments(args, profiler: FeatureQuantifier):
    AlnRunner = {
        "bwa": BwaMemRunner,
        "minimap2": Minimap2Runner,
    }.get(args.aligner)

    if AlnRunner is None:
        raise ValueError(f"Aligner `{args.aligner}` is not supported.")

    aln_runner = AlnRunner(
        args.cpus_for_alignment,
        args.reference,
        sample_id=os.path.basename(args.out_prefix),
    )

    for i, (input_type, *reads) in enumerate(args.input_data):

        logger.info("Running %s alignment: %s", input_type, ",".join(reads))

        sam_suffix = f".{input_type}.{i}"

        debug_samfile, samfile = None, None
        if args.keep_alignment_file:
            samfile = args.keep_alignment_file.replace(".sam", "")
            samfile = f"{samfile}{sam_suffix}.sam"
        if profiler.debug:
            debug_samfile = f"{profiler.out_prefix}{sam_suffix}.filtered.sam"

        proc, call = aln_runner.run(reads, single_end_reads=input_type == "single", alignment_file=samfile)

        # if proc.returncode != 0:
        #     logger.error("Encountered problems aligning.")
        #     logger.error("Aligner call was:")
        #     logger.error("%s", call)
        #     logger.error("Shutting down.")
        #     sys.exit(1)

        # pylint: disable=W0718
        try:

            profiler.count_alignments(
                proc.stdout,
                aln_format="sam",
                min_identity=args.min_identity,
                min_seqlen=args.min_seqlen,
                sam_prefix=f".{input_type}.{i}",
                debug_samfile=debug_samfile,
            )

        except Exception as err:
            if isinstance(err, ValueError) and str(err).strip() == "file does not contain alignment data":
                # pylint: disable=W1203
                logger.error(f"Failed to align. Is `{args.aligner}` installed and on the path?")
                logger.error("Aligner call was:")
                logger.error("%s", call)
                raise ValueError from err

            logger.error("Encountered problems digesting the alignment stream:")
            logger.error("%s", err)
            logger.error("Aligner call was:")
            logger.error("%s", call)
            logger.error("Shutting down.")
            raise ValueError from err


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
    logger.info("Command: %s %s", os.path.basename(sys.argv[0]), " ".join(sys.argv[1:]))

    kwargs = {"debug": args.debug}
    annotation_db = args.annotation_db
    db_input = None
    if args.run_mode == RunMode.GENE:
        Quantifier = GeneQuantifier
    else:
        Quantifier, kwargs["run_mode"] = RegionQuantifier, args.run_mode
        db_args = {}
        if args.run_mode == RunMode.DOMAIN:
            annotation_db = SmallDatabaseImporter(
                single_category="feature", db_format=args.db_format,
            )
            db_input = args.annotation_db
        elif args.run_mode == RunMode.SMALL_GENOME:
            annotation_db = SmallGenomeDatabaseImporter()
            try:
                db_input, db_args["input_data2"] = args.annotation_db.split(",")
            except ValueError as exc:
                raise ValueError("Require two input files.") from exc
        else:
            kwargs["in_memory"] = args.db_in_memory

        if db_input:
            annotation_db.build_database(
                db_input,
                **db_args,
            )
        logger.info("Finished loading database.")

    profiler = Quantifier(
        db=annotation_db,
        out_prefix=args.out_prefix,
        distribution_mode=args.distribution_mode,
        strand_specific=args.strand_specific,
        paired_end_count=args.paired_end_count,
        calculate_coverage=args.with_coverage,
        **kwargs,
    )

    if args.input_type == "fastq":

        stream_alignments(args, profiler)

    else:

        input_file = args.bam if args.input_type == "bam" else args.sam
        debug_samfile = None
        if profiler.debug:
            debug_samfile = f"{profiler.out_prefix}.{args.input_type}.filtered.sam"

        profiler.count_alignments(
            sys.stdin if input_file == "-" else input_file,
            aln_format=args.input_type,
            min_identity=args.min_identity,
            min_seqlen=args.min_seqlen,
            external_readcounts=args.import_readcounts,
            unmarked_orphans=args.unmarked_orphans,
            debug_samfile=debug_samfile,
        )

    profiler.finalise(
        restrict_reports=args.restrict_metrics,
        report_category=True,
        report_unannotated=args.run_mode.report_unannotated,
        dump_counters=args.debug,
        in_memory=args.db_in_memory,
        gene_group_db=args.gene_group_db,
    )


if __name__ == "__main__":
    main()

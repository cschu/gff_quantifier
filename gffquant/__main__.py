# pylint: disable=C0103,C0301

""" module docstring """

import logging
import os
import pathlib
import sys

# pylint: disable=W0611
from .db.db_import import SmallDatabaseImporter
from .handle_args import handle_args
from .ui.validation import check_input_reads
from .profilers import GeneQuantifier, RegionQuantifier
from .runners.alignment_runner import BwaMemRunner, Minimap2Runner
from . import __version__


logger = logging.getLogger(__name__)


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
    logger.info("Command: %s %s", os.path.basename(sys.argv[0]), " ".join(sys.argv[1:]))

    kwargs = {}
    annotation_db = args.annotation_db
    if args.mode in ("gene", "genes"):
        Quantifier = GeneQuantifier
    else:
        Quantifier, kwargs["reference_type"] = RegionQuantifier, args.mode
        if args.mode == "domain":
            annotation_db = SmallDatabaseImporter(
                logger, args.annotation_db, single_category="feature", sep=args.db_separator,
            )
            logger.info("Finished loading database.")

    input_data = check_input_reads(
        fwd_reads=args.reads1, rev_reads=args.reads2,
        single_reads=args.singles, orphan_reads=args.orphans,
    )

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
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
            proc, call = aln_runner.run(reads, single_end_reads=input_type == "single", alignment_file=args.keep_alignment_file)

            # if proc.returncode != 0:
            #     logger.error("Encountered problems aligning.")
            #     logger.error("Aligner call was:")
            #     logger.error("%s", call)
            #     logger.error("Shutting down.")
            #     sys.exit(1)

            # pylint: disable=W0718
            try:

                profiler.count_alignments(
                    proc.stdout, aln_format="sam", min_identity=args.min_identity, min_seqlen=args.min_seqlen,
                )

            except Exception as err:
                if isinstance(err, ValueError) and str(err).strip() == "file does not contain alignment data":
                    # pylint: disable=W1203
                    logger.error(f"Failed to align. Is `{args.aligner}` installed and on the path?")
                    logger.error("Aligner call was:")
                    logger.error("%s", call)
                    sys.exit(1)

                logger.error("Encountered problems digesting the alignment stream:")
                logger.error("%s", err)
                logger.error("Aligner call was:")
                logger.error("%s", call)
                logger.error("Shutting down.")
                sys.exit(1)





    else:
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

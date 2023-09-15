# pylint: disable=C0103,R0914

""" module docstring """

import gzip
import json
import logging
import os
import time

from abc import ABC, abstractmethod
from collections import Counter

from ..alignment import AlignmentGroup, AlignmentProcessor, SamFlags
from ..annotation import GeneCountAnnotator, RegionCountAnnotator, CountWriter
from ..counters import CountManager
from ..db.annotation_db import AnnotationDatabaseManager

from .. import __tool__


logger = logging.getLogger(__name__)


class FeatureQuantifier(ABC):
    # pylint: disable=R0902,R0913
    TRUE_AMBIG_MODES = ("dist1", "1overN")

    def __init__(
        self,
        db=None,
        out_prefix=__tool__,
        ambig_mode="unique_only",
        reference_type="genome",
        strand_specific=False,
        paired_end_count=1,
    ):
        self.aln_counter = Counter()
        self.db = db
        self.adm = None
        self.do_overlap_detection = reference_type in ("genome", "domain")
        self.reference_type = reference_type
        self.count_manager = CountManager(
            distribution_mode=ambig_mode,
            region_counts=reference_type in ("genome", "domain"),
            strand_specific=strand_specific and reference_type not in ("genome", "domain"),
            paired_end_count=paired_end_count,
        )
        self.out_prefix = out_prefix
        self.ambig_mode = ambig_mode
        self.bamfile = None
        self.alp = None
        self.reference_manager = {}
        self.strand_specific = strand_specific

    def allow_ambiguous_alignments(self):
        """All distribution modes support ambiguous alignments,
        with the exception of `unique_only."""
        return self.ambig_mode != "unique_only"

    def treat_ambiguous_as_unique(self):
        """In 'non-true' ambig modes (all1, primary_only),
        ambiguous alignments are treated as unique alignments."""
        return self.ambig_mode not in FeatureQuantifier.TRUE_AMBIG_MODES

    def require_ambig_bookkeeping(self):
        """'Ambiguous bookkeeping' is a special treatment for true ambig modes.
        It requires all alignments of a read to be processed together."""
        return (
            self.allow_ambiguous_alignments() and not self.treat_ambiguous_as_unique()
        )

    def process_alignments_sameref(self, ref, alignments, aln_count=1):
        """
        Process a group of alignments against the same reference sequence.

        The weird output format is due to compatibility reasons with ambiguous alignments
        when one of the TRUE_AMBIG_MODES is chosen.

        input:
          - reference clear text identifier (required for overlap detection)
          - list of alignments (reference short id, start, end, is_reverse)

        output:
          - generator! with the following items:
            in overlap mode:
              {reference_id:
                  [(begin overlap, end overlap,
				    is_reverse, begin covered region, end covered region)
               for each overlapped region]}
                  + alignment count: 1 if the alignment overlaps any region of interest, else 0
                  + unalignment count: 1 - alignment count
            else:
                  {reference_id:(-1, -1, is_reverse, -1, -1)
                  + alignment count: 1
				  (each reference is already a feature, i.e. everything is aligned)
                  + unalignment count: 0
        """

        for rid, start, end, rev_strand in alignments:
            if self.do_overlap_detection:

                hits = {
                    (sstart, send, rev_strand, None, None)
                    for (sstart, send)
                    in self.adm.get_overlaps(
                        ref, start, end,
                        domain_mode=self.reference_type == "domain",
                    )
                }

                # if the alignment overlaps multiple features, each one gets a count
                aln_count = int(bool(hits)) * aln_count

            else:
                hits = {(None, None, rev_strand, None, None)}
                # aln_count = 1 / aln_count

            yield ({rid: hits}, aln_count, 0 if aln_count else 1)

    def process_counters(
        self,
        restrict_reports=None,
        report_category=True,
        report_unannotated=True,
        dump_counters=True,
    ):
        if self.adm is None:
            self.adm = AnnotationDatabaseManager.from_db(self.db)

        if dump_counters:
            self.count_manager.dump_raw_counters(self.out_prefix, self.reference_manager)

        report_scaling_factors = restrict_reports is None or "scaled" in restrict_reports

        annotator_type = (GeneCountAnnotator, RegionCountAnnotator)[self.do_overlap_detection]
        count_annotator = annotator_type(self.strand_specific, report_scaling_factors=report_scaling_factors)
        count_annotator.annotate(self.reference_manager, self.adm, self.count_manager)

        count_writer = CountWriter(
            self.out_prefix,
            has_ambig_counts=self.count_manager.has_ambig_counts(),
            strand_specific=self.strand_specific,
            restrict_reports=restrict_reports,
            report_category=report_category,
            total_readcount=self.aln_counter["read_count"],
            filtered_readcount=self.aln_counter["filtered_read_count"],
        )

        unannotated_reads = self.count_manager.get_unannotated_reads()
        unannotated_reads += self.aln_counter["unannotated_ambig"]

        count_writer.write_feature_counts(
            self.adm,
            count_annotator,
            (None, unannotated_reads)[report_unannotated],
        )

        count_writer.write_gene_counts(
            count_annotator.gene_counts,
            count_annotator.scaling_factors["total_gene_uniq"],
            count_annotator.scaling_factors["total_gene_ambi"]
        )

        self.adm.clear_caches()

    def register_reference(self, rid, aln_reader):
        known_ref = self.reference_manager.get(rid)
        new_ref = aln_reader.get_reference(rid)

        if known_ref is None or new_ref == known_ref:
            self.reference_manager[rid] = new_ref
        else:
            raise ValueError(f"Reference clash {rid} points to old_ref={known_ref} and new_ref={new_ref}.")

        return new_ref[0]

    @abstractmethod
    def process_alignment_group(self, aln_group, aln_reader):
        ...

    def process_alignments(self, aln_reader, min_identity=None, min_seqlen=None, unmarked_orphans=False):
        # pylint: disable=R0914
        t0 = time.time()

        aln_stream = aln_reader.get_alignments(
            min_identity=min_identity,
            min_seqlen=min_seqlen,
            allow_multiple=self.allow_ambiguous_alignments(),
            allow_unique=True,
            filter_flags=SamFlags.SUPPLEMENTARY_ALIGNMENT,
        )

        self.count_manager.toggle_single_read_handling(unmarked_orphans)

        aln_count = 0
        read_count = 0
        current_aln_group = None
        for aln_count, aln in enumerate(aln_stream, start=1):
            if self.ambig_mode == "primary_only" and not aln.is_primary():
                continue
            if self.ambig_mode in ("uniq_only", "unique_only") and not aln.is_unique():
                continue

            if current_aln_group is None or current_aln_group.qname != aln.qname:
                if current_aln_group is not None:
                    self.process_alignment_group(current_aln_group, aln_reader)
                current_aln_group = AlignmentGroup()
                read_count += 1

                if read_count and read_count % 100000 == 0:
                    logger.info("Processed %s reads.", read_count)

            current_aln_group.add_alignment(aln)

        if current_aln_group is not None:
            self.process_alignment_group(current_aln_group, aln_reader)

        if aln_count == 0:
            logger.warning("No alignments present in stream.")

        t1 = time.time()
        logger.info("Processed %s reads (%s alignments) in %s.", read_count, aln_count, f"{t1 - t0:.3f}s")

        return aln_count, read_count, 0, None

    @staticmethod
    def get_readcount(internal_readcounts, external_readcounts, verbose=True):
        # pylint: disable=W0703
        # need to figure out what exceptions to catch...
        read_count = internal_readcounts
        if os.path.isfile(external_readcounts):
            try:
                with open(external_readcounts, encoding="UTF-8") as read_counts_in:
                    read_count = json.load(read_counts_in)["n_reads"]
                if verbose:
                    logger.info("Found pre-filter readcounts (%s).", read_count)
            except Exception as err:
                print(f"Error accessing readcounts: {err}")
                logger.warning(
                    "Could not access pre-filter readcounts. Using post-filter readcounts (%s).\n"
                    "This should result in an alignment-rate of 100%%.",
                    read_count
                )
        else:
            read_count = int(external_readcounts)

        return read_count

    # pylint: disable=W0613
    def count_alignments(
        self,
        aln_stream,
        aln_format="sam",
        min_identity=None,
        min_seqlen=None,
        external_readcounts=None,
        unmarked_orphans=False,
    ):
        aln_reader = AlignmentProcessor(aln_stream, aln_format)

        aln_count, _, unannotated_ambig, _ = self.process_alignments(
            aln_reader,
            min_identity=min_identity,
            min_seqlen=min_seqlen,
            unmarked_orphans=unmarked_orphans,
        )

        full_readcount, read_count, filtered_readcount = aln_reader.read_counter

        if external_readcounts is not None:
            full_readcount = external_readcounts
            # FeatureQuantifier.get_readcount(full_read_count, external_readcounts)

        self.aln_counter.update(
            {
                "aln_count": aln_count,
                "read_count": read_count,
                "unannotated_ambig": unannotated_ambig,
                "full_read_count": full_readcount,
                "filtered_read_count": filtered_readcount,
            }
        )

        self.aln_counter.update(aln_reader.get_alignment_stats_dict())

    def finalise(
        self,
        restrict_reports=None,
        report_category=False,
        report_unannotated=False,
        dump_counters=False,
    ):

        with gzip.open(f"{self.out_prefix}.aln_stats.txt.gz", "wt") as aln_stats_out:
            print(
                AlignmentProcessor.get_alignment_stats_str(
                    [
                        v
                        for k, v in self.aln_counter.items()
                        if k.startswith("pysam_") and not k.endswith("total")
                    ],
                    table=True,
                ),
                file=aln_stats_out
            )

        if self.aln_counter.get("aln_count"):
            self.process_counters(
                restrict_reports=restrict_reports,
                report_category=report_category,
                report_unannotated=report_unannotated,
                dump_counters=dump_counters,
            )

            # for metric, value in self.aln_counter.items():
            #     logger.info("%s: %s", metric, value)

            for metric, value in (
                ("Input reads", "full_read_count"),
                ("Aligned reads", "read_count"),
                ("Alignments", "pysam_total"),
                ("Reads passing filters", "filtered_read_count"),
                ("Alignments passing filters", "pysam_passed"),
                ("  - Discarded due to seqid", "pysam_seqid_filt"),
                ("  - Discarded due to length", "pysam_len_filt"),
                # ("Unannotated multimappers", "unannotated_ambig"),
            ):
                logger.info("%s: %s", metric, self.aln_counter.get(value))

            """
            aln_count: 6469971
            [2023-08-30 16:51:23,259] read_count: 2498819
            [2023-08-30 16:51:23,259] unannotated_ambig: 0
            [2023-08-30 16:51:23,259] full_read_count: 2574299
            [2023-08-30 16:51:23,260] filtered_read_count: 2241860
            [2023-08-30 16:51:23,260] pysam_total: 11875545
            [2023-08-30 16:51:23,260] pysam_passed: 6469971
            [2023-08-30 16:51:23,260] pysam_seqid_filt: 1683770
            [2023-08-30 16:51:23,260] pysam_len_filt: 3721804
            """

            logger.info(
                "Alignment rate: %s%%, Filter pass rate: %s%%",
                round(self.aln_counter["read_count"] / self.aln_counter["full_read_count"], 3) * 100,
                round(self.aln_counter["filtered_read_count"] / self.aln_counter["full_read_count"], 3) * 100,
            )

        logger.info("Finished.")

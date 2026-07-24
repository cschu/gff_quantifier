# pylint: disable=C0103,R0914,W1203,R0902,R0917

""" module docstring """

import json
import logging
import os
import time

from abc import ABC
from collections import Counter

import numpy as np

from .panda_coverage_profiler import PandaCoverageProfiler
from ..alignment import AlignmentGroup, AlignmentProcessor, ReferenceHit, SamFlags
from ..annotation import GeneCountAnnotator, RegionCountAnnotator, CountWriter
from ..counters import AlignmentCounter
from ..counters.count_matrix import CountMatrix
from ..db.annotation_db import AnnotationDatabaseManager
from ..util.fileio import get_open_function

from .. import __tool__, DistributionMode, RunMode


logger = logging.getLogger(__name__)


class FeatureQuantifier(ABC):
    """
        Three groups of alignments:
        1. Ambiguous alignments ('multimappers')
        - single-end
        - input orphans
        - newly orphaned (one mate didn't align)
        - complete pairs
        2. Properly aligning pairs, which have to be
        + unique (each mate aligns once)
        + concordant (mates align to the same reference)
        3. Other uniques. Here we could have uniquely aligning
        - single-end reads
        - input orphan reads
        - newly orphaned reads, e.g. if only one mate of an input pair aligns
        - pairs, with mates aligning to different references
            -> i.e. reads could be derived from either kind of library
        - these could also be multimappers in all1-mode!
        """
    # pylint: disable=R0902,R0913

    def __init__(
        self,
        db=None,
        out_prefix=__tool__,
        distribution_mode=DistributionMode.ONE_OVER_N,
        run_mode=RunMode.GENE,
        strand_specific=False,
        paired_end_count=1,
        calculate_coverage=False,
        external_genecounts=None,
        debug=False,
    ):
        self.aln_counter = Counter()
        self.db = db
        self.adm = None
        self.run_mode = run_mode
        self.counter = AlignmentCounter(
            distribution_mode=distribution_mode,
            strand_specific=strand_specific,
            paired_end_count=paired_end_count,
            initial_cols=6 if external_genecounts else 2,
        )
        self.out_prefix = out_prefix
        self.distribution_mode = distribution_mode
        self.reference_manager = {}
        self.strand_specific = strand_specific
        self.debug = debug
        self.panda_cv = PandaCoverageProfiler(dump_dataframes=self.debug) if calculate_coverage else None

        self.external_counts = bool(external_genecounts)
        self.total_reads = None
        self.filtered_reads = None
        if external_genecounts:
            self._import_counts(external_genecounts)

    def _import_counts(self, fn):
        with get_open_function(fn)(fn, "rt") as _in:
            try:
                _ = next(_in)
            except StopIteration as exc:
                raise ValueError(f"Counts file is empty: {fn}") from exc

            gene_id, *counts = next(_in).strip().split("\t")
            ncols = len(counts)
            self.total_reads = float(next(_in).split("\t")[0])
            self.filtered_reads = float(next(_in).split("\t")[0])

            for i, row in enumerate(_in, start=1,):
                gene_id, *counts = row.strip().split("\t")
                # counts = np.array(tuple(map(float, counts)), dtype=np.float64)
                self.reference_manager[i] = (gene_id, 1)
                self.counter[i] = np.array(tuple(map(float, counts)), dtype=CountMatrix.NUMPY_DTYPE)

        return self.counter.counts.colsums()

    def check_hits(self, ref, aln):
        """ Check if an alignment hits a region of interest on a reference sequence.
            - in overlap modes, this performs an overlap check against the annotation database
            - in gene mode, every alignment is a hit
            - Returns
              - `has_target`, indicating if the reference sequence contains one or more potential hits
              - list of hits (alignment can overlap multiple regions)
        """
        if self.run_mode.overlap_required:
            overlaps = self.adm.get_overlaps(
                ref, aln.start, aln.end,
                domain_mode=self.run_mode == RunMode.DOMAIN,
                with_coverage=self.run_mode == RunMode.SMALL_GENOME,
            )

            # interrogate overlap stream header if overlaps were found
            has_target = next(overlaps).has_target

            hits = [
                ReferenceHit(
                    rid=aln.rid,
                    start=ovl_target.start,
                    end=ovl_target.end,
                    rev_strand=aln.is_reverse(),
                    cov_start=ovl_target.cov_start,
                    cov_end=ovl_target.cov_end,
                    has_annotation=ovl_target.has_annotation(),
                )
                for ovl_target in overlaps
            ]

        else:

            # if no overlap-check was performed, each hit is a hit on target (e.g. aligning against gene catalogues)
            has_target, hits = True, [ReferenceHit(rid=aln.rid, rev_strand=aln.is_reverse())]

        return has_target, hits

    def process_counters(
        self,
        restrict_reports=None,
        report_category=True,
        report_unannotated=True,
        dump_counters=True,
        in_memory=True,
        gene_group_db=False,
    ):
        if self.adm is None:
            self.adm = AnnotationDatabaseManager.from_db(self.db, in_memory=in_memory)

        if dump_counters and not self.external_counts:
            self.counter.dump(self.out_prefix, self.reference_manager,)

        report_scaling_factors = restrict_reports is None or "scaled" in restrict_reports

        if self.run_mode.overlap_required and not self.external_counts:
            Annotator = RegionCountAnnotator
        else:
            Annotator = GeneCountAnnotator

        count_annotator = Annotator(self.strand_specific, report_scaling_factors=report_scaling_factors)

        if self.external_counts:
            logger.info("TOTAL_GENE_COUNTS = %s (IMPORTED)", self.counter.counts.colsums())
        else:
            total_gene_counts = self.counter.generate_gene_count_matrix(self.reference_manager)
            logger.info("TOTAL_GENE_COUNTS = %s", total_gene_counts)
            self.total_reads = self.aln_counter["read_count"]
            self.filtered_reads = self.aln_counter["filtered_read_count"]

        count_writer = CountWriter(
            self.out_prefix,
            has_ambig_counts=self.counter.has_ambig_counts(),
            strand_specific=self.strand_specific,
            restrict_reports=restrict_reports,
            report_category=report_category,
            total_readcount=self.total_reads,
            filtered_readcount=self.filtered_reads,
        )

        if not self.external_counts or self.debug:
            count_writer.write_gene_counts(
                self.counter,
                self.reference_manager,
                self.total_reads,
                self.filtered_reads,
            )

        ggroups = tuple(
            (self.reference_manager.get(key[0] if isinstance(key, tuple) else key))[0]  # .split(".")[0]
            for key, _ in self.counter
        )

        if dump_counters or self.debug:
            self.counter.counts.dump(prefix=self.out_prefix, labels=ggroups)

        self.counter.group_gene_count_matrix(self.reference_manager)
        unannotated_reads = self.counter.counts[CountMatrix.NO_ANNOTATION] + self.aln_counter["unannotated_ambig"]

        if dump_counters or self.debug:
            self.counter.counts.dump(prefix=self.out_prefix, state="ggroup")

        functional_counts, category_sums = count_annotator.annotate_gene_counts(
            self.reference_manager,
            self.adm,
            self.counter,
            gene_group_db=gene_group_db,
        )

        logger.info("FC-index: %s", str(list(functional_counts.index.keys())[:10]))
        logger.info("FC-counts: %s", str(functional_counts.counts[0:10, :]))

        categories = self.adm.get_categories()
        for category, category_sum in zip(categories, category_sums):
            features = tuple(self.adm.get_features(category.id))
            rows = tuple(
                key[0] == category.id
                for key, _ in functional_counts
            )

            cat_counts = CountMatrix.from_count_matrix(functional_counts, rows=rows)

            logger.info("PROCESSING CATEGORY=%s", category.name)
            count_writer.write_category(
                category.id,
                category.name,
                category_sum,
                cat_counts,
                features,
                unannotated_reads=(None, unannotated_reads)[report_unannotated],
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

    def process_alignments(
        self,
        aln_reader: AlignmentProcessor,
        sam_prefix="",
        min_identity=None,
        min_seqlen=None,
        unmarked_orphans=False,
        debug_samfile=None,
        # panda: PandaProfiler=None,
    ):
        # pylint: disable=R0914
        t0 = time.time()

        aln_stream = aln_reader.get_alignments(
            min_identity=min_identity,
            min_seqlen=min_seqlen,
            filter_flags=SamFlags.SUPPLEMENTARY_ALIGNMENT,
            filtered_sam=debug_samfile,
        )

        self.counter.toggle_single_read_handling(unmarked_orphans)
        ac = self.aln_counter

        read_count = 0
        current_aln_group = None

        for ac["aln_count"], aln in enumerate(aln_stream, start=1):

            has_target, aln.hits = self.check_hits(aln.refname, aln)
            # do we care about alignments to unannotated regions?
            # for quantification purposes i'd say, we don't, but maybe
            # this needs to be a user-defined option
            if not has_target:
                continue

            self.reference_manager.setdefault(aln.rid, (aln.refname, aln.reflength))

            ac["alignments_on_target"] += 1

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

        if ac["aln_count"] == 0:
            logger.warning("No alignments present in stream.")

        t1 = time.time()
        logger.info("Processed %s reads (%s alignments) in %s.", read_count, ac["aln_count"], f"{t1 - t0:.3f}s")

        return ac["aln_count"], read_count, 0, None

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
        sam_prefix="",
        debug_samfile=None,
    ):
        aln_reader = AlignmentProcessor(aln_stream, aln_format)

        self.process_alignments(
            aln_reader,
            min_identity=min_identity,
            min_seqlen=min_seqlen,
            unmarked_orphans=unmarked_orphans,
            sam_prefix=sam_prefix,
            debug_samfile=debug_samfile,
        )

        full_readcount, read_count, filtered_readcount = aln_reader.read_counter

        if external_readcounts is not None:
            full_readcount = external_readcounts

        self.aln_counter.update(
            {
                "read_count": read_count,
                "unannotated_ambig": 0,  # did this come from region analysis, ie. ambig reads that don't hit an annotated region?
                "full_read_count": full_readcount,
                "filtered_read_count": filtered_readcount,
            }
        )

        self.aln_counter.update(aln_reader.get_alignment_stats_dict())

    def report_alignments(self):
        with open(f"{self.out_prefix}.aln_stats.json", "wt", encoding="UTF-8",) as aln_stats_out:
            json.dump(self.aln_counter, aln_stats_out)

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

        if self.aln_counter["full_read_count"]:
            nreads = self.aln_counter["read_count"]
            nallreads = self.aln_counter["full_read_count"]
            nfiltered = self.aln_counter["filtered_read_count"]
            alignment_rate = round(nreads / nallreads, 3) * 100
            filter_pass_rate = round(nfiltered / nallreads, 3) * 100
        else:
            alignment_rate, filter_pass_rate = None, None

        logger.info(
            "Alignment rate: %s%%, Filter pass rate: %s%%",
            alignment_rate, filter_pass_rate,
        )

    def finalise(
        self,
        restrict_reports=None,
        report_category=False,
        report_unannotated=False,
        dump_counters=False,
        in_memory=True,
        gene_group_db=False,
    ):

        if self.aln_counter.get("aln_count") or self.external_counts:
            if self.adm is None:
                self.adm = AnnotationDatabaseManager.from_db(self.db, in_memory=in_memory)

            # self.write_coverage()

            self.process_counters(
                restrict_reports=restrict_reports,
                report_category=report_category,
                report_unannotated=report_unannotated,
                dump_counters=dump_counters,
                in_memory=in_memory,
                gene_group_db=gene_group_db,
            )

            self.adm.clear_caches()

        logger.info("Finished.")

    def evaluate_alignment_group(self, aln_group, is_ambig_alignment):
        ...

    def process_alignment_group(self, aln_group, aln_reader):
        """ Checks an alignment group for hits and passes the hit list
        to the counters.
        """
        is_ambiguous_group = aln_group.is_ambiguous()
        is_ambig_alignment = self.distribution_mode.require_ambig_tracking and is_ambiguous_group

        keep_group = True

        if is_ambiguous_group and not self.distribution_mode.allow_ambiguous:
            # if the alignment group has ambiguous alignments, but we don't allow them: get out
            keep_group = False

        elif not self.distribution_mode.allow_secondary:
            # if we only allow primary alignments, clear the secondaries
            _ = [lst.clear() for lst in aln_group.secondaries]

            keep_group = bool(aln_group.primaries)

        if keep_group:
            # # hit counts are mode-dependent:
            # # ambiguous alignments in true ambig modes require the number of alternative locations
            # # (we handle first/second PE-alignments separately, hence a two-element list)
            # # otherwise the hit count is always 1
            # ambig_hit_counts = [0, 0]
            # all_hits = []

            # # process each alignment of the group individually
            # for aln in aln_group.get_alignments():
            #     current_ref = self.register_reference(aln.rid, aln_reader)
            #     # check for region hits (in overlap modes, otherwise alignments always hit)
            #     has_target, hits = self.check_hits(current_ref, aln)
            #     if has_target and hits:
            #         if is_ambig_alignment:
            #             # update the hit counts
            #             ambig_hit_counts[aln.is_second()] += 1
            #         all_hits.append((aln, hits))

            # # pass the hit lists to the counters
            # count_stream = (
            #     (aln_hits, ambig_hit_counts[aln.is_second()] if is_ambig_alignment else 1)
            #     for aln, aln_hits in all_hits
            # )
            count_stream = tuple(
                aln_group.get_all_hits(
                    as_ambiguous=self.distribution_mode.require_ambig_tracking
                )
            )

            contributed_counts = self.counter.update(
                count_stream,
                ambiguous_counts=is_ambiguous_group,
                pair=aln_group.is_paired(),
                pe_library=aln_group.pe_library,
            )

            ambig_hit_counts = aln_group.ambig_hit_counts

            msg = f"{ambig_hit_counts=} {contributed_counts=})"

        else:

            msg = f"--> DROPPED distribution_mode={self.distribution_mode.alias}"

        logger.debug(
            f"aln_group {aln_group.qname} "
            f"(ambig={is_ambig_alignment} size={aln_group.n_align()} "
            f"{msg}"
        )

        return keep_group

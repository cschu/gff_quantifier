# pylint: disable=C0103,R0914,W1203

""" module docstring """

import gzip
import json
import logging
import os
import time

from abc import ABC
from collections import Counter
from dataclasses import dataclass

import pandas as pd

from ..alignment import AlignmentGroup, AlignmentProcessor, SamFlags
from ..annotation import GeneCountAnnotator, RegionCountAnnotator, CountWriter
from ..counters import CountManager
from ..db.annotation_db import AnnotationDatabaseManager

from .. import __tool__, DistributionMode, RunMode


logger = logging.getLogger(__name__)


@dataclass
class ReferenceHit:
    rid: int = None
    start: int = None
    end: int = None
    rev_strand: bool = None
    cov_start: int = None
    cov_end: int = None
    has_annotation: bool = None
    n_aln: int = None
    is_ambiguous: bool = None
    library_mod: int = None

    def __hash__(self):
        return hash(tuple(self.__dict__.values()))

    def __eq__(self, other):
        return all(
            item[0][1] == item[1][1]
            for item in zip(
                sorted(self.__dict__.items()),
                sorted(other.__dict__.items())
            )
        )

    def __str__(self):
        return "\t".join(map(str, self.__dict__.values()))

    def __repr__(self):
        return str(self)


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
    ):
        self.aln_counter = Counter()
        self.db = db
        self.adm = None
        self.run_mode = run_mode
        self.count_manager = CountManager(
            distribution_mode=distribution_mode,
            region_counts=run_mode.overlap_required,
            strand_specific=strand_specific and not run_mode.overlap_required,
            paired_end_count=paired_end_count,
        )
        self.out_prefix = out_prefix
        self.distribution_mode = distribution_mode
        self.reference_manager = {}
        self.strand_specific = strand_specific
        self.coverage_counter = {}

    def update_coverage(self, aln_hits):
        for hits, n_aln in aln_hits:
            for hit in hits:
                self.coverage_counter.setdefault(hit.is_ambiguous, {}).setdefault((hit.rid, hit.start, hit.end), Counter()).update({p: 1 / n_aln for p in range(hit.cov_start, hit.cov_end)})
    def _calc_coverage(self):
        for key in sorted(set(self.coverage_counter.get(True, {})).union(self.coverage_counter.get(False, {}))):
            uniq_cov, ambig_cov = self.coverage_counter.get(True, {}).get(key, Counter()), self.coverage_counter.get(False, {}).get(key, Counter())
            length = key[2] - key[1] + 1
            len_both = len(set(uniq_cov).union(ambig_cov))
            yield {
                "rid": key[0],
                "start": key[1],
                "end": key[2],
                "length": length,
                "uniq_depth": sum(uniq_cov) / length,
                "uniq_depth_covered": (sum(uniq_cov) / len(uniq_cov)) if uniq_cov else 0.0,
                "uniq_horizontal": len(uniq_cov) / length,
                "combined_depth": (sum(uniq_cov) + sum(ambig_cov)) / length,
                "combined_depth_covered": ((sum(uniq_cov) + sum(ambig_cov)) / len_both) if len_both else 0.0,
                "combined_horizontal": len_both / length,
            }
    def write_coverage(self):
        df = pd.DataFrame(self._calc_coverage())
        
        categories = {
            cat.id: cat
            for cat in self.adm.get_categories()
        }
# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# OOOOOOOOO              OOOOOOOO
# OOOOOOOOO              OOOOOOOO
# OOOOOOOOO              OOOOOOOO
# OOOOOOOOO              OOOOOOOO
# OOOOOOOOO              OOOOOOOO
# OOOOOOOOO              OOOOOOOO
# OOOOOOOOO              OOOOOOOO
# OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        annotated_cols = []
        for rid, start, end in zip(df["rid"], df["start"], df["end"]):
            ref, reflen = self.reference_manager.get(rid)
            for annseq in self.adm.get_db_sequence(ref, start=start, end=end):
                if annseq.annotation_str is not None: #and start == annseq.start and annseq.end == end:
                    d = {"refid": rid, "start": start, "end": end, "refname": annseq.featureid}
                    d.update({cat.name: None for cat in categories.values()})
                    annotated_cols.append(d)
                    for item in annseq.annotation_str.split(";"):
                        catid, features = item.split("=")
                        d[categories.get(int(catid)).name] = [int(feat) for feat in features.split(",")]

        df2 = pd.DataFrame.from_records(annotated_cols)
        coverage_columns = ["uniq_depth", "uniq_depth_covered", "uniq_horizontal", "combined_depth", "combined_depth_covered", "combined_horizontal"]

        for category in categories.values():
            features = pd.DataFrame.from_records(
                {"fid": feat.id, "feature": feat.name}
                for feat in self.adm.get_features(category=category.id)
                # {feat.id: feat.name for feat in self.adm.get_features(category=category.id)}
            )

            cat_grouped = pd.merge(
                df2[["refid", "start", "end", "refname", category.name]],
                df,
                left_index=False, right_index=False,
                left_on=("refid", "start", "end",),
                right_on=("rid", "start", "end",),
            ) \
                .dropna(axis=0, subset=[category.name,], how="any") \
                .explode(category.name, ignore_index=True)[[category.name,] + coverage_columns] \
                .groupby(category.name, as_index=False)
            # df3 = pd.merge(df1[["refid","start","end","refname","PFAMs"]], df2, left_index=False, right_index=False, left_on=("refid", "start", "end"), right_on=("rid","start","end"))

            # >>> df1.explode("PFAMs").groupby(by="PFAMs")[["start", "PFAMs"]].mean("start")
            coverage_df = cat_grouped[[category.name, "uniq_horizontal", "combined_horizontal",]].mean(numeric_only=True)
            depth_df = cat_grouped[[category.name, "uniq_depth", "uniq_depth_covered", "combined_depth", "combined_depth_covered",]].sum(numeric_only=True)            
            
            out_df = pd.merge(
                features,
                # cat_grouped.mean(),
                pd.merge(coverage_df, depth_df, on=(category.name,), left_index=False, right_index=False),
                left_index=False,
                right_index=False,
                left_on=("fid",),
                right_on=(category.name,),
            ).drop([category.name, "fid"], axis=1)

            new_order = ["feature", "uniq_depth", "uniq_depth_covered", "uniq_horizontal", "combined_depth", "combined_depth_covered", "combined_horizontal",]
            out_df[new_order].to_csv(
                f"{self.out_prefix}.{category.name}.coverage.txt", sep="\t", index=False, float_format="%.5f"
            )


            # cat_grouped.mean().to_csv(f"{self.out_prefix}.{category.name}.coverage.txt", sep="\t", index=False, float_format="%.5f")


        df.to_csv(self.out_prefix + ".all.coverage.txt", index=False, sep="\t", na_rep="NA")
        df2.to_csv(self.out_prefix + ".all.coverage_annotation.txt", index=False, sep="\t", na_rep="NA")

        self.adm.dump(self.out_prefix + ".db")
    

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

            # has_target, *_ = next(overlaps)
            # interrogate overlap stream header if overlaps were found
            has_target = next(overlaps).has_target

            hits = [
                ReferenceHit(
                    rid=aln.rid,
                    # start=start,
                    start=ovl_target.start,
                    # end=end,
                    end = ovl_target.end,
                    rev_strand=aln.is_reverse(),
                    cov_start=ovl_target.cov_start,
                    cov_end=ovl_target.cov_end,
                    # has_annotation=(dbseq is not None and dbseq.annotation_str),
                    has_annotation=ovl_target.has_annotation(),                    
                )
                # for _, start, end, dbseq in overlaps
                for ovl_target in overlaps
            ]

        else:

            has_target, hits = True, [ReferenceHit(rid=aln.rid, rev_strand=aln.is_reverse())]

        return has_target, hits

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

        Annotator = (GeneCountAnnotator, RegionCountAnnotator)[self.run_mode.overlap_required]
        count_annotator = Annotator(self.strand_specific, report_scaling_factors=report_scaling_factors)
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

    def process_alignments(self, aln_reader, min_identity=None, min_seqlen=None, unmarked_orphans=False):
        # pylint: disable=R0914
        t0 = time.time()

        aln_stream = aln_reader.get_alignments(
            min_identity=min_identity,
            min_seqlen=min_seqlen,
            filter_flags=SamFlags.SUPPLEMENTARY_ALIGNMENT,
            filtered_sam=f"{self.out_prefix}.filtered.sam",
        )

        self.count_manager.toggle_single_read_handling(unmarked_orphans)
        ac = self.aln_counter

        read_count = 0
        current_aln_group = None
        for ac["aln_count"], aln in enumerate(aln_stream, start=1):

            has_target, aln.hits = self.check_hits(aln.refname, aln)
            # do we care about alignments to unannotated regions?
            # for quantification purposes i'd say, we don't, but maybe
            # this needs to be a user-defined option
            if not has_target:
                continue

            self.reference_manager.setdefault(aln.rid, (aln.refname, aln.reflength))
            
            ac["alignments_on_target"] += 1

            if current_aln_group is None or current_aln_group.qname != aln.qname:
                if current_aln_group is not None:
                    for hit in self.process_alignment_group(current_aln_group, aln_reader):
                        yield hit
                current_aln_group = AlignmentGroup()
                read_count += 1

                if read_count and read_count % 100000 == 0:
                    logger.info("Processed %s reads.", read_count)

            current_aln_group.add_alignment(aln)

        if current_aln_group is not None:
            for hit in self.process_alignment_group(current_aln_group, aln_reader):
                yield hit

        if ac["aln_count"] == 0:
            logger.warning("No alignments present in stream.")

        t1 = time.time()
        logger.info("Processed %s reads (%s alignments) in %s.", read_count, ac["aln_count"], f"{t1 - t0:.3f}s")



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

        hits = self.process_alignments(
            aln_reader,
            min_identity=min_identity,
            min_seqlen=min_seqlen,
            unmarked_orphans=unmarked_orphans,
        )

        # pd.DataFrame(hits).to_csv(self.out_prefix + ".hits.tsv", sep="\t", index=False)
        raw_df = pd.DataFrame(hits)
        raw_df.to_csv(self.out_prefix + ".hits.tsv", sep="\t", index=False)
        raw_df["contrib"] = 1 / raw_df["n_aln"] / raw_df["library_mod"]

        keep_columns = ["rid", "start", "end", "contrib"]
        contrib_sums_uniq = raw_df[raw_df["is_ambiguous"] == False][keep_columns].groupby(by=["rid", "start", "end"], as_index=False).sum()
        contrib_sums_combined = raw_df[keep_columns].groupby(by=["rid", "start", "end"], as_index=False).sum()
        raw_df = pd.merge(
            contrib_sums_uniq,
            contrib_sums_combined,
            on=("rid", "start", "end"),
            left_index=False, right_index=False,
            how="outer",
        ).rename({"contrib_x": "uniq_raw", "contrib_y": "combined_raw"}, axis=1).fillna(0)
        raw_df["uniq_lnorm"] = raw_df["uniq_raw"] / (raw_df["end"] - raw_df["start"] + 1)
        raw_df["combined_lnorm"] = raw_df["combined_raw"] / (raw_df["end"] - raw_df["start"] + 1)
        raw_df.to_csv(self.out_prefix + ".raw_lnorm.tsv", sep="\t", index=False)
        
        self.write_coverage()

        full_readcount, read_count, filtered_readcount = aln_reader.read_counter

        if external_readcounts is not None:
            full_readcount = external_readcounts            

        self.aln_counter.update(
            {                
                "read_count": read_count,
                "unannotated_ambig": 0,
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

            logger.info(
                "Alignment rate: %s%%, Filter pass rate: %s%%",
                round(self.aln_counter["read_count"] / self.aln_counter["full_read_count"], 3) * 100,
                round(self.aln_counter["filtered_read_count"] / self.aln_counter["full_read_count"], 3) * 100,
            )

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
            self.update_coverage(count_stream)


            contributed_counts = self.count_manager.update_counts(
                count_stream,
                ambiguous_counts=is_ambiguous_group,
                pair=aln_group.is_paired(),
                pe_library=aln_group.pe_library,
            )
            ambig_hit_counts = aln_group.ambig_hit_counts
            msg = f"{ambig_hit_counts=} {contributed_counts=})"            

            for hits, n_aln in count_stream:
                for hit in hits:
                    hit.n_aln = n_aln
                    yield hit

        else:

            msg = f"--> DROPPED distribution_mode={self.distribution_mode.alias}"

        logger.debug(
            f"aln_group {aln_group.qname} "
            f"(ambig={is_ambig_alignment} size={aln_group.n_align()} "
            f"{msg}"
        )

        # return keep_group

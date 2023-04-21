# pylint: disable=C0103,W1514,R0913

""" module docstring """

import gzip
import logging

import numpy as np


logger = logging.getLogger(__name__)


class CountWriter:
    COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled", "rpkm"]

    def __init__(
        self,
        prefix,
        aln_count,
        has_ambig_counts=False,
        strand_specific=False,
        restrict_reports=None,
        report_category=True,
        report_unannotated=True,
    ):
        self.out_prefix = prefix
        self.aln_count = aln_count
        self.has_ambig_counts = has_ambig_counts
        self.strand_specific = strand_specific
        self.publish_reports = [
            item for item in CountWriter.COUNT_HEADER_ELEMENTS
            if restrict_reports is None or item in restrict_reports
        ]
        if report_category:
            self.publish_reports.append("category")
        if report_unannotated:
            self.publish_reports.append("unannotated")

    def get_header(self):
        reports = [
            report
            for report in self.publish_reports
            if report in CountWriter.COUNT_HEADER_ELEMENTS
        ]
        header = []
        header += (f"uniq_{element}" for element in reports)
        if self.has_ambig_counts:
            header += (
                f"combined_{element}" for element in reports
            )
        if self.strand_specific:
            for strand in ("ss", "as"):
                header += (
                    f"uniq_{element}_{strand}"
                    for element in reports
                )
                if self.has_ambig_counts:
                    header += (
                        f"combined_{element}_{strand}"
                        for element in reports
                    )
        return header

    def compile_output_row(self, counts, scaling_factor=1, ambig_scaling_factor=1):

        def compile_block(raw, lnorm, scaling_factors):
            return (raw, lnorm,) + tuple(lnorm * factor for factor in scaling_factors)

        p, row = 0, []
        rpkm_factor = 1e9 / self.aln_count
        # unique counts
        row += compile_block(*counts[p:p + 2], (scaling_factor, rpkm_factor,))
        p += 2
        # ambiguous counts
        if self.has_ambig_counts:
            row += compile_block(*counts[p:p + 2], (ambig_scaling_factor, rpkm_factor,))
            p += 2
        # sense-strand unique
        if self.strand_specific:
            row += compile_block(*counts[p:p + 2], (scaling_factor, rpkm_factor,))
            p += 2
            # sense-strand ambiguous
            if self.has_ambig_counts:
                row += compile_block(*counts[p:p + 2], (ambig_scaling_factor, rpkm_factor,))
                p += 2
            # antisense-strand unique
            row += compile_block(*counts[p:p + 2], (scaling_factor, rpkm_factor,))
            p += 2
            if self.has_ambig_counts:
                row += compile_block(*counts[p:p + 2], (ambig_scaling_factor, rpkm_factor,))

        out_row = []
        n = len(CountWriter.COUNT_HEADER_ELEMENTS)
        for col, item in enumerate(row):
            if CountWriter.COUNT_HEADER_ELEMENTS[col % n] in self.publish_reports:
                out_row.append(item)

        return out_row

    def write_feature_counts(self, db, unannotated_reads, featcounts):
        for category_id, counts in sorted(featcounts.items()):
            scaling_factor, ambig_scaling_factor = featcounts.scaling_factors[
                category_id
            ]
            category = db.query_category(category_id).name
            if "scaled" in self.publish_reports:
                logger.info(
                    "SCALING FACTORS %s %s %s",
                    category, scaling_factor, ambig_scaling_factor
                )
            with gzip.open(f"{self.out_prefix}.{category}.txt.gz", "wt") as feat_out:
                print("feature", *self.get_header(), sep="\t", file=feat_out)
                if "unannotated" in self.publish_reports:
                    print("unannotated", unannotated_reads, sep="\t", file=feat_out)

                if "category" in self.publish_reports:
                    cat_counts = counts.get(f"cat:::{category_id}")
                    if cat_counts is not None:
                        out_row = self.compile_output_row(
                            cat_counts,
                            scaling_factor=featcounts.scaling_factors["total_uniq"],
                            ambig_scaling_factor=featcounts.scaling_factors["total_ambi"],
                        )
                        print(
                            "category",
                            *(f"{c:.5f}" for c in out_row),
                            flush=True,
                            sep="\t",
                            file=feat_out,
                        )

                for feature_id, f_counts in sorted(counts.items()):
                    if feature_id.startswith("cat:::"):
                        continue
                    feature = db.query_feature(feature_id).name
                    out_row = self.compile_output_row(
                        f_counts,
                        scaling_factor=scaling_factor,
                        ambig_scaling_factor=ambig_scaling_factor,
                    )
                    print(
                        feature,
                        *(f"{c:.5f}" for c in out_row),
                        flush=True,
                        sep="\t",
                        file=feat_out,
                    )

    def write_gene_counts(self, gene_counts, uniq_scaling_factor, ambig_scaling_factor):
        if "scaled" in self.publish_reports:
            logger.info("SCALING_FACTORS %s %s", uniq_scaling_factor, ambig_scaling_factor)
        with gzip.open(f"{self.out_prefix}.gene_counts.txt.gz", "wt") as gene_out:
            print("gene", *self.get_header(), sep="\t", file=gene_out, flush=True)

            for gene, g_counts in sorted(gene_counts.items()):
                out_row = self.compile_output_row(
                    g_counts,
                    scaling_factor=uniq_scaling_factor,
                    ambig_scaling_factor=ambig_scaling_factor
                )
                print(gene, *(f"{c:.5f}" for c in out_row), flush=True, sep="\t", file=gene_out)

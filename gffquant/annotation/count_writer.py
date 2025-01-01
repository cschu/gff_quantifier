# pylint: disable=C0103,W1514,R0913,R0917,R0914

""" module docstring """

import gzip
import logging
import sys

import numpy as np

from ..counters import AlignmentCounter


logger = logging.getLogger(__name__)


class CountWriter:
    COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled", "rpkm"]

    def __init__(
        self,
        prefix,
        has_ambig_counts=False,
        strand_specific=False,
        restrict_reports=None,
        report_category=True,
        total_readcount=None,
        filtered_readcount=None,
    ):
        self.out_prefix = prefix
        self.has_ambig_counts = has_ambig_counts
        self.strand_specific = strand_specific
        self.publish_reports = [
            item for item in CountWriter.COUNT_HEADER_ELEMENTS
            if restrict_reports is None or item in restrict_reports
        ]
        if report_category:
            self.publish_reports.append("category")
        if total_readcount:
            self.publish_reports.append("total_readcount")
        if filtered_readcount:
            self.publish_reports.append("filtered_readcount")

        self.total_readcount = total_readcount
        self.filtered_readcount = filtered_readcount

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
        rpkm_factor = 1e9 / self.filtered_readcount

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

    @staticmethod
    def write_row(header, data, stream=sys.stdout):
        print(header, *(f"{c:.5f}" for c in data), flush=True, sep="\t", file=stream)

    def write_category(
        self,
        category_id,
        category_name,
        category_sum,
        counts,
        feature_names,
        unannotated_reads=None,
        report_unseen=True,
    ):
        with gzip.open(f"{self.out_prefix}.{category_name}.txt.gz", "wt") as feat_out:
            header = self.get_header()
            print("feature", *header, sep="\t", file=feat_out)

            if unannotated_reads is not None:
                print("unannotated", unannotated_reads, sep="\t", file=feat_out)

            if "total_readcount" in self.publish_reports:
                CountWriter.write_row(
                    "total_reads",
                    np.zeros(len(header)) + self.total_readcount,
                    stream=feat_out,
                )

            if "filtered_readcount" in self.publish_reports:
                CountWriter.write_row(
                    "filtered_reads",
                    np.zeros(len(header)) + self.filtered_readcount,
                    stream=feat_out,
                )

            if "category" in self.publish_reports:
                # cat_counts = counts[0]
                cat_counts = category_sum
                logger.info("CAT %s: %s", category_name, str(cat_counts))
                if cat_counts is not None:
                    CountWriter.write_row("category", category_sum, stream=feat_out)

            # for item in counts:
            #     if not isinstance(item[0], tuple):
            #         logger.info("ITEM: %s", str(item))
            #         raise TypeError(f"Weird key: {str(item)}")
            #     (cid, fid), fcounts = item
            #     if (report_unseen or fcounts.sum()) and cid == category_id:
            #         CountWriter.write_row(feature_names[fid], fcounts, stream=feat_out,)

            for (cid, fid), fcounts in counts:
                if (report_unseen or fcounts.sum()) and cid == category_id:
                    CountWriter.write_row(feature_names[fid], fcounts, stream=feat_out,)

    def write_gene_counts(
        self,
        gene_counts: AlignmentCounter,
        refmgr,
        gene_group_db=False,
    ):
        with gzip.open(f"{self.out_prefix}.gene_counts.txt.gz", "wt") as gene_out:
            print("gene", *self.get_header(), sep="\t", file=gene_out, flush=True)

            ref_stream = (
                (
                    refmgr.get(rid[0] if isinstance(rid, tuple) else rid)[0],
                    rid,
                )
                for rid, _ in gene_counts
            )

            for ref, rid in sorted(ref_stream):
                counts = gene_counts[rid]
                if gene_group_db:
                    ref_tokens = ref.split(".")
                    gene_id, _ = ".".join(ref_tokens[:-1]), ref_tokens[-1]
                else:
                    gene_id = ref

                CountWriter.write_row(gene_id, counts, stream=gene_out,)

# pylint: disable=C0103,W1514,R0913,R0917

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

    def write_category(self, category, counts, index, names, unique_sf, ambig_sf, unannotated_reads=None, report_unseen=True):
        # category, c_counts, c_index, c_names, u_sf, a_sf
        if "scaled" in self.publish_reports:
            logger.info(
                "SCALING FACTORS %s %s %s",
                category, unique_sf, ambig_sf,
            )
        with gzip.open(f"{self.out_prefix}.{category}.txt.gz", "wt") as feat_out:
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
                # cat_counts = counts.get(f"cat:::{category_id}")
                cat_counts = counts[0]  # np.array((counts[0][0], counts[0][2], counts[0][1], counts[0][3]))
                logger.info("CAT %s: %s", category, str(cat_counts))
                if cat_counts is not None:
                    cat_row = self.compile_output_row(
                        cat_counts,
                        # scaling_factor=featcounts.scaling_factors["total_uniq"],
                        # ambig_scaling_factor=featcounts.scaling_factors["total_ambi"],
                        scaling_factor=unique_sf,
                        ambig_scaling_factor=ambig_sf,
                    )
                    CountWriter.write_row("category", cat_row, stream=feat_out)

            for fid, i in index.items():
                f_counts = counts[i]  # np.array((counts[i][0], counts[i][2], counts[i][1], counts[i][3]))  #counts[fid]
                if report_unseen or f_counts.sum():
                    out_row = self.compile_output_row(
                        f_counts,
                        scaling_factor=unique_sf,
                        ambig_scaling_factor=ambig_sf,
                    )
                    CountWriter.write_row(names[fid], out_row, stream=feat_out)


    # pylint: disable=R0914
    def write_feature_counts(self, db, featcounts, unannotated_reads=None, report_unseen=True):
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
                    cat_counts = counts.get(f"cat:::{category_id}")
                    logger.info("CAT %s: %s", category_id, str(cat_counts))
                    if cat_counts is not None:
                        cat_row = self.compile_output_row(
                            cat_counts,
                            scaling_factor=featcounts.scaling_factors["total_uniq"],
                            ambig_scaling_factor=featcounts.scaling_factors["total_ambi"],
                        )
                        CountWriter.write_row("category", cat_row, stream=feat_out)

                for feature in db.get_features(category_id):
                    f_counts = counts.get(str(feature.id), np.zeros(len(header)))
                    if report_unseen or f_counts.sum():
                        out_row = self.compile_output_row(
                            f_counts,
                            scaling_factor=scaling_factor,
                            ambig_scaling_factor=ambig_scaling_factor,
                        )
                        CountWriter.write_row(feature.name, out_row, stream=feat_out)

    def write_gene_counts(
        self,
        gene_counts: AlignmentCounter,
        refmgr,
        uniq_scaling_factor,
        ambig_scaling_factor,
        gene_group_db=False
    ):
        if "scaled" in self.publish_reports:
            logger.info("SCALING_FACTORS %s %s", uniq_scaling_factor, ambig_scaling_factor)
        with gzip.open(f"{self.out_prefix}.gene_counts.txt.gz", "wt") as gene_out:
            print("gene", *self.get_header(), sep="\t", file=gene_out, flush=True)

            ref_stream = (
                (
                    refmgr.get(rid[0] if isinstance(rid, tuple) else rid)[0],
                    rid,
                )
                for rid in gene_counts.get_all_regions()
            )

            for ref, rid in sorted(ref_stream):
                counts = gene_counts.get_counts(rid)
                if gene_group_db:
                    ref_tokens = ref.split(".")
                    gene_id, _ = ".".join(ref_tokens[:-1]), ref_tokens[-1]
                else:
                    gene_id = ref

                out_row = self.compile_output_row(
                    counts,
                    scaling_factor=uniq_scaling_factor,
                    ambig_scaling_factor=ambig_scaling_factor,
                )

                CountWriter.write_row(gene_id, out_row, stream=gene_out,)

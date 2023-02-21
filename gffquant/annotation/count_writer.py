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
        reports = self.publish_reports
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

    # pylint: disable=R0914
    def write_coverage(self, db, coverage_counts):
        uniq_cov, combined_cov = {}, {}
        uniq_depth, combined_depth = {}, {}
        uniq_depth_raw, combined_depth_raw = {}, {}
        counts = (
            uniq_cov, combined_cov,
            uniq_depth_raw, combined_depth_raw,
            uniq_depth, combined_depth
        )

        all_categories, all_features = set(), set()
        for cov_data in coverage_counts.values():
            n_features = sum(len(features) for _, features in cov_data["annotation"])
            for category, features in cov_data["annotation"]:
                all_categories.add(category)
                all_features.update(features)
                for feature in features:
                    # coverage = number of positions that are covered by at least 1 alignment
                    if cov_data["uniq_coverage"].any():
                        uniq_cov.setdefault(category, {}).setdefault(feature, []).append(
                            (cov_data["uniq_coverage"][cov_data["uniq_coverage"] > 0]).mean()
                        )
                    if cov_data["combined_coverage"].any():
                        combined_cov.setdefault(category, {}).setdefault(feature, []).append(
                            (cov_data["combined_coverage"][cov_data["combined_coverage"] > 0]).mean()
                        )
                    # raw depth = sum of the read depths over all positions;
                    # in case of a region annotated with more than one feature,
                    # raw depth and depth are divided by number of features
                    uniq_depth_raw.setdefault(category, {}).setdefault(feature, []).append(
                        cov_data["uniq_coverage"].sum() / n_features
                    )
                    combined_depth_raw.setdefault(category, {}).setdefault(feature, []).append(
                        cov_data["combined_coverage"].sum() / n_features
                    )
                    # depth = average read depth over all positions
                    uniq_depth.setdefault(category, {}).setdefault(feature, []).append(
                        cov_data["uniq_coverage"].mean() / n_features
                    )
                    combined_depth.setdefault(category, {}).setdefault(feature, []).append(
                        cov_data["combined_coverage"].mean() / n_features
                    )

        for category_id in sorted(all_categories):
            category = db.query_category(category_id).name
            with gzip.open(f"{self.out_prefix}.{category}.coverage.txt.gz", "wt") as feat_out:
                print(
                    "category",
                    "feature",
                    "coverage_unique",
                    "coverage_combined",
                    "depth_raw_unique",
                    "depth_raw_combined",
                    "depth_unique",
                    "depth_combined",

                    sep="\t",
                    file=feat_out,
                )
                for feature_id in sorted(all_features):
                    feature = db.query_feature(feature_id).name
                    uc, cc, urd, crd, ud, cd = (
                        np.sum(d.get(category_id, {}).get(feature_id, (0.0,)))
                        for d in counts
                    )

                    print(
                        category, feature,
                        uc, cc, urd, crd, ud, cd,
                        sep="\t", file=feat_out,
                    )


# pylint: disable=W0105
"""
domain_cov = dict()
        for pos in set(ref_coverage).union(ambig_ref_coverage):
            # features_at_position = pos_feature_dict.get(pos, set())
            for domtype in pos_feature_dict.get(pos, set()):
                domain_cov.setdefault(domtype, dict()).update(
                    {
                        "depth_uniq": list(),
                        "depth_ambig": list(),
                        "cov_uniq": list(),
                        "cov_ambig": list(),
                    }
                )

                uniq_cov = ref_coverage.get(pos, Counter())
                ambig_cov = ambig_ref_coverage.get(pos, Counter())

                if uniq_cov:
                    domain_cov[domtype]["depth_uniq"].append(
                        sum(uniq_cov.values()) / len(uniq_cov.values())
                    )
                    domain_cov[domtype]["cov_uniq"].append(
                        sum(1 for v in uniq_cov.values() if v) / len(uniq_cov.values())
                    )
                if ambig_cov:
                    domain_cov[domtype]["depth_ambig"].append(
                        sum(ambig_cov.values()) / len(ambig_cov.values())
                    )
                    domain_cov[domtype]["cov_ambig"].append(
                        sum(1 for v in ambig_cov.values() if v) / len(ambig_cov.values())
                    )

        with gzip.open(self.out_prefix + ".covsum.txt.gz", "wt") as cov_out:
            print("#domain", "depth_unique", "depth_combined", "coverage_unique", "coverage_combined", sep="\t", file=cov_out,)
            for domtype, counts in sorted(domain_cov.items()):

                depth_uniq, depth_ambig = [
                    c for c in counts.get("depth_uniq", list()) if c is not None
                ], [c for c in counts.get("depth_ambig", list()) if c is not None]

                depth_uniq_ = (sum(depth_uniq) / len(depth_uniq)) if depth_uniq else 0
                print(domtype, "UNIQ", depth_uniq, sum(depth_uniq), len(depth_uniq), "=", depth_uniq_,)

                depth_ambig_ = (
                    (sum(depth_ambig) / len(depth_ambig)) if depth_ambig else 0
                )
                print(domtype, "AMBIG", depth_ambig, sum(depth_ambig), len(depth_ambig), "=", depth_ambig_,)

                if depth_ambig_ < depth_uniq_:
                    raise ValueError(
                        f"{domtype}: depth_uniq_ + depth_ambig_ {depth_ambig_:.5f} is smaller than uniq depth {depth_uniq_:.5f}."
                    )

                cov_uniq, cov_ambig = [
                    c for c in counts.get("cov_uniq", list()) if c is not None
                ], [c for c in counts.get("cov_ambig", list()) if c is not None]

                cov_uniq_ = (sum(cov_uniq) / len(cov_uniq)) if cov_uniq else 0
                print(
                    domtype,
                    "UNIQ",
                    cov_uniq,
                    sum(cov_uniq),
                    len(cov_uniq),
                    "=",
                    cov_uniq_,
                )

                cov_ambig_ = (sum(cov_ambig) / len(cov_ambig)) if cov_ambig else 0
                print(
                    domtype,
                    "AMBIG",
                    cov_ambig,
                    sum(cov_ambig),
                    len(cov_ambig),
                    "=",
                    cov_ambig_,
                )

                if cov_ambig_ < cov_uniq_:
                    raise ValueError(
                        f"{domtype}: cov_uniq_ + cov_ambig_ {cov_ambig_:.5f} is smaller than uniq cov {cov_uniq_:.5f}."
                    )

                # print(domtype, f"{depth_uniq_:.5f}", f"{depth_ambig_:.5f}", sep="\t", file=cov_out)
                print(
                    domtype,
                    *(
                        f"{val:.5f}"
                        for val in (depth_uniq_, depth_ambig_, cov_uniq_, cov_ambig_)
                    ),
                    sep="\t",
                    file=cov_out,
                )"""

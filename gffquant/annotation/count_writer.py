# pylint: disable=C0103,W1514

""" module docstring """

import gzip

import numpy as np


class CountWriter:
    COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled"]

    def __init__(self, prefix, has_ambig_counts=False, strand_specific=False):
        self.out_prefix = prefix
        self.has_ambig_counts = has_ambig_counts
        self.strand_specific = strand_specific

    def get_header(self):
        header = []
        header += (f"uniq_{element}" for element in CountWriter.COUNT_HEADER_ELEMENTS)
        if self.has_ambig_counts:
            header += (
                f"combined_{element}" for element in CountWriter.COUNT_HEADER_ELEMENTS
            )
        if self.strand_specific:
            for strand in ("ss", "as"):
                header += (
                    f"uniq_{element}_{strand}"
                    for element in CountWriter.COUNT_HEADER_ELEMENTS
                )
                if self.has_ambig_counts:
                    header += (
                        f"combined_{element}_{strand}"
                        for element in CountWriter.COUNT_HEADER_ELEMENTS
                    )
        return header

    def compile_output_row(self, counts, scaling_factor=1, ambig_scaling_factor=1):
        p, row = 0, []
        # unique counts
        row += tuple(counts[p:p + 2])
        row += (row[-1] * scaling_factor,)
        p += 2
        # ambiguous counts
        if self.has_ambig_counts:
            row += tuple(counts[p:p + 2])
            row += (row[-1] * ambig_scaling_factor,)
            p += 2
        # sense-strand unique
        if self.strand_specific:
            row += tuple(counts[p:p + 2])
            row += (row[-1] * scaling_factor,)
            p += 2
            # sense-strand ambiguous
            if self.has_ambig_counts:
                row += tuple(counts[p:p + 2])
                row += (row[-1] * ambig_scaling_factor,)
                p += 2
            # antisense-strand unique
            row += tuple(counts[p:p + 2])
            row += (row[-1] * scaling_factor,)
            p += 2
            if self.has_ambig_counts:
                row += tuple(counts[p:p + 2])
                row += (row[-1] * ambig_scaling_factor,)
        return row  # + [scaling_factor, ambig_scaling_factor]

    def write_feature_counts(self, db, unannotated_reads, featcounts):
        for category_id, counts in sorted(featcounts.items()):
            scaling_factor, ambig_scaling_factor = featcounts.scaling_factors[
                category_id
            ]
            category = db.query_category(category_id).name
            print("SCALING FACTORS", category, scaling_factor, ambig_scaling_factor)
            with gzip.open(f"{self.out_prefix}.{category}.txt.gz", "wt") as feat_out:
                print("feature", *self.get_header(), sep="\t", file=feat_out)
                print("unannotated", unannotated_reads, sep="\t", file=feat_out)
                for feature_id, f_counts in sorted(counts.items()):
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
        print("SCALING_FACTORS", uniq_scaling_factor, ambig_scaling_factor)
        with gzip.open(f"{self.out_prefix}.gene_counts.txt.gz", "wt") as gene_out:
            print("gene", *self.get_header(), sep="\t", file=gene_out, flush=True)

            for gene, g_counts in sorted(gene_counts.items()):
                print(gene, g_counts)
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

        all_categories, all_features = set(), set()
        for cov_data in coverage_counts.values():
            n_features = sum(len(features) for _, features in cov_data["annotation"])
            for category, features in cov_data["annotation"]:
                all_categories.add(category)
                all_features.update(features)
                for feature in features:
                    uniq_depth.setdefault(category, {}).setdefault(feature, []).append(
                        cov_data["uniq_coverage"].mean() / n_features
                    )
                    combined_depth.setdefault(category, {}).setdefault(feature, []).append(
                        cov_data["combined_coverage"].mean() / n_features
                    )
                    uniq_cov.setdefault(category, {}).setdefault(feature, []).append(
                        (cov_data["uniq_coverage"] > 0).mean()
                    )
                    combined_cov.setdefault(category, {}).setdefault(feature, []).append(
                        (cov_data["combined_coverage"] > 0).mean()
                    )

        for category_id in sorted(all_categories):
            category = db.query_category(category_id).name
            with gzip.open(f"{self.out_prefix}.{category}.coverage.txt.gz", "wt") as feat_out:
                print(
                    "category",
                    "feature",
                    "coverage_unique",
                    "coverage_combined",
                    "depth_unique",
                    "depth_combined",
                    sep="\t",
                    file=feat_out,
                )
                for feature_id in sorted(all_features):
                    feature = db.query_feature(feature_id).name
                    uc, cc, ud, cd = [
                        d.get(category_id, {}).get(feature_id)
                        for d in (uniq_cov, combined_cov, uniq_depth, combined_depth)
                    ]

                    print(
                        category,
                        feature,
                        np.mean(uc) if uc is not None else "NA",
                        np.mean(cc) if cc is not None else "NA",
                        np.mean(ud) if ud is not None else "NA",
                        np.mean(cd) if cd is not None else "NA",
                        sep="\t",
                        file=feat_out,
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

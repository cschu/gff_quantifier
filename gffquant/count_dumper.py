# pylint: disable=C0103,W1514

""" module docstring """


class CountDumper:
    COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled"]

    def __init__(self, prefix, has_ambig_counts=False, strand_specific=False):
        self.out_prefix = prefix
        self.has_ambig_counts = has_ambig_counts
        self.strand_specific = strand_specific

    def get_header(self):
        header = []
        header += (f"uniq_{element}" for element in CountDumper.COUNT_HEADER_ELEMENTS)
        if self.has_ambig_counts:
            header += (
                f"combined_{element}" for element in CountDumper.COUNT_HEADER_ELEMENTS
            )
        if self.strand_specific:
            for strand in ("ss", "as"):
                header += (
                    f"uniq_{element}_{strand}"
                    for element in CountDumper.COUNT_HEADER_ELEMENTS
                )
                if self.has_ambig_counts:
                    header += (
                        f"combined_{element}_{strand}"
                        for element in CountDumper.COUNT_HEADER_ELEMENTS
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
        return row

    def dump_feature_counts(self, unannotated_reads, featcounts):
        with open(f"{self.out_prefix}.feature_counts.txt", "w") as feat_out:
            print(
                "category",
                "feature",
                *self.get_header(),
                sep="\t",
                file=feat_out,
                flush=True,
            )
            print(
                "unannotated",
                "",
                unannotated_reads,
                sep="\t",
                file=feat_out,
                flush=True,
            )
            for category, counts in sorted(featcounts.items()):
                scaling_factor, ambig_scaling_factor = featcounts.scaling_factors[
                    category
                ]
                for feature, f_counts in sorted(counts.items()):
                    out_row = self.compile_output_row(
                        f_counts,
                        scaling_factor=scaling_factor,
                        ambig_scaling_factor=ambig_scaling_factor,
                    )
                    print(
                        category,
                        feature,
                        *(f"{c:.5f}" for c in out_row),
                        flush=True,
                        sep="\t",
                        file=feat_out,
                    )

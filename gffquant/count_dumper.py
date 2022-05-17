# pylint: disable=C0103,W1514

""" module docstring """

import gzip


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
        return row + [scaling_factor, ambig_scaling_factor]

    def dump_feature_counts(self, unannotated_reads, featcounts):
        for category, counts in sorted(featcounts.items()):
            scaling_factor, ambig_scaling_factor = featcounts.scaling_factors[
                category
            ]
            with gzip.open(f"{self.out_prefix}.{category}.txt.gz", "wt") as feat_out:
                print("feature", *self.get_header(), sep="\t", file=feat_out)
                print("unannotated", unannotated_reads, sep="\t", file=feat_out)
                for feature, f_counts in sorted(counts.items()):
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

    def dump_gene_counts(self, gene_counts, uniq_scaling_factor, ambig_scaling_factor):
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

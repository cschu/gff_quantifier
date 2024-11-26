""" module docstring """

import pickle

from collections import Counter

import pandas as pd

from .panda_profiler import PandaProfiler


# # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
# # OOOOOOOOO              OOOOOOOO
# # OOOOOOOOO              OOOOOOOO
# # OOOOOOOOO              OOOOOOOO
# # OOOOOOOOO              OOOOOOOO
# # OOOOOOOOO              OOOOOOOO
# # OOOOOOOOO              OOOOOOOO
# # OOOOOOOOO              OOOOOOOO
# # OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO


class PandaCoverageProfiler(PandaProfiler):
    COLUMNS = [
        "uniq_depth",
        "uniq_depth_covered",
        "uniq_horizontal",
        "combined_depth",
        "combined_depth_covered",
        "combined_horizontal",
    ]

    def __init__(self, dump_dataframes=False):
        PandaProfiler.__init__(self, with_overlap=True)
        self._coverage_data = [{}, {}]
        self.main_df = None
        self.dump_dataframes = dump_dataframes
        self.mode = "coverage"

    def update_coverage(self, aln_hits):
        for hits, n_aln in aln_hits:
            for hit in hits:
                cov_key = hit.cov_start, hit.cov_end
                hit_key = hit.rid, hit.start, hit.end
                increment = (1 / hit.library_mod) / n_aln
                self._coverage_data[hit.is_ambiguous] \
                    .setdefault(hit_key, Counter())[cov_key] += increment
        # for hits, n_aln in aln_hits:
        #     for hit in hits:
        #         self._coverage_data[hit.is_ambiguous] \
        #             .setdefault((hit.rid, hit.start, hit.end), Counter()) \
        #             .update({p: (1 / hit.library_mod / n_aln) for p in range(hit.cov_start, hit.cov_end + 1)})

    def _calc_coverage(self):
        for key in sorted(
            set(self._coverage_data[False]).union(self._coverage_data[True])
        ):
            uniq_cov, ambig_cov = Counter(), Counter()
            for (cstart, cend), counts in self._coverage_data[False].get(key, Counter()).items():
                uniq_cov.update({p: counts for p in range(cstart, cend + 1)})
            for (cstart, cend), counts in self._coverage_data[True].get(key, Counter()).items():
                ambig_cov.update({p: counts for p in range(cstart, cend + 1)})

            length = key[2] - key[1] + 1
            len_both = len(set(uniq_cov).union(ambig_cov))

            sum_uniq_cov, sum_ambig_cov = sum(uniq_cov.values()), sum(ambig_cov.values())

            yield {
                "rid": key[0],
                "start": key[1],
                "end": key[2],
                "length": length,
                "uniq_depth": sum_uniq_cov / length,
                "uniq_depth_covered": (sum_uniq_cov / len(uniq_cov)) if uniq_cov else 0.0,
                "uniq_horizontal": len(uniq_cov) / length,
                "combined_depth": (sum_uniq_cov + sum_ambig_cov) / length,
                "combined_depth_covered": ((sum_uniq_cov + sum_ambig_cov) / len_both) if len_both else 0.0,
                "combined_horizontal": len_both / length,
            }

        # for key in sorted(
        #     # set(self._coverage_data.get(True, {})) \
        #     #     .union(self._coverage_data.get(False, {}))
        #     set(self._coverage_data[False]).union(self._coverage_data[True])
        # ):
        #     uniq_cov = self._coverage_data[False].get(key, Counter())
        #     ambig_cov = self._coverage_data[True].get(key, Counter())
        #     length = key[2] - key[1] + 1
        #     len_both = len(set(uniq_cov).union(ambig_cov))

        #     sum_uniq_cov, sum_ambig_cov = sum(uniq_cov.values()), sum(ambig_cov.values())

        #     yield {
        #         "rid": key[0],
        #         "start": key[1],
        #         "end": key[2],
        #         "length": length,
        #         "uniq_depth": sum_uniq_cov / length,
        #         "uniq_depth_covered": (sum_uniq_cov / len(uniq_cov)) if uniq_cov else 0.0,
        #         "uniq_horizontal": len(uniq_cov) / length,
        #         "combined_depth": (sum_uniq_cov + sum_ambig_cov) / length,
        #         "combined_depth_covered": ((sum_uniq_cov + sum_ambig_cov) / len_both) if len_both else 0.0,
        #         "combined_horizontal": len_both / length,
        #     }

    def dump_coverage(self, read_data_provider, out_prefix):
        if self.dump_dataframes:
            with open(f"{read_data_provider.out_prefix}.coverage.dat", "wb") as _out:
                pickle.dump(self._coverage_data, _out)
        self.main_df = pd.DataFrame(self._calc_coverage())
        self.main_df.to_csv(out_prefix + ".main_before_annotation.txt", index=False, sep="\t", na_rep="NA")

        self._annotate_records(
            self.get_gene_coords(),
            read_data_provider.reference_manager,
            read_data_provider.adm,
        )

        categories = {cat.id: cat for cat in read_data_provider.adm.get_categories()}

        gene_category_map = self._get_gene_category_map(categories, read_data_provider)
        if self.dump_dataframes:
            self.main_df.to_csv(out_prefix + ".all.coverage.txt", index=False, sep="\t", na_rep="NA")
            gene_category_map.to_csv(out_prefix + ".all.coverage_annotation.txt", index=False, sep="\t", na_rep="NA")
            read_data_provider.adm.dump(out_prefix + ".db")

        for category in categories.values():
            features = pd.DataFrame.from_records(
                {"fid": feat.id, "feature": feat.name}
                for feat in read_data_provider.adm.get_features(category=category.id)
            )

            columns = [category.name,] + PandaCoverageProfiler.COLUMNS

            cat_grouped = self._annotate_category_counts(
                self.main_df,
                gene_category_map,
                category.name,
            ) \
                .explode(category.name, ignore_index=True)[columns] \
                .groupby(category.name, as_index=False)

            coverage_df = cat_grouped[
                [category.name, "uniq_horizontal", "combined_horizontal",]
            ] \
                .mean(numeric_only=True)
            depth_df = cat_grouped[
                [
                    category.name,
                    "uniq_depth",
                    "uniq_depth_covered",
                    "combined_depth",
                    "combined_depth_covered",
                ]
            ] \
                .sum(numeric_only=True)

            out_df = pd.merge(
                features,
                pd.merge(coverage_df, depth_df, on=(category.name,), left_index=False, right_index=False),
                left_index=False,
                right_index=False,
                left_on=("fid",),
                right_on=(category.name,),
            ) \
                .drop([category.name, "fid"], axis=1) \
                .sort_values(by=["feature",])

            new_order = ["feature"] + PandaCoverageProfiler.COLUMNS
            out_df[new_order].to_csv(
                f"{out_prefix}.{category.name}.coverage.txt", sep="\t", index=False, float_format="%.5f"
            )

        gene_columns = ["gene"] + PandaCoverageProfiler.COLUMNS
        self.main_df[gene_columns] \
            .sort_values(by=["gene",]) \
            .to_csv(out_prefix + ".genes.coverage.txt", index=False, sep="\t", na_rep="NA", float_format="%.5f")

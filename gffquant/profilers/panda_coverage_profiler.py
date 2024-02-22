from collections import Counter

import pandas as pd

from .panda_profiler import PandaProfiler


class PandaCoverageProfiler(PandaProfiler):
    def __init__(self):
        PandaProfiler.__init__(self)
        self._coverage_data = {}
        self.main_df = None

def update_coverage(self, aln_hits):
    for hits, n_aln in aln_hits:
        for hit in hits:
            self._coverage_data.setdefault(hit.is_ambiguous, {}) \
                .setdefault((hit.rid, hit.start, hit.end), Counter()) \
                .update({p: 1 / n_aln for p in range(hit.cov_start, hit.cov_end)})

def _calc_coverage(self):
        for key in sorted(
            set(self._coverage_data.get(True, {})) \
                .union(self._coverage_data.get(False, {}))
        ):
            uniq_cov = self._coverage_data.get(True, {}).get(key, Counter())
            ambig_cov = self.coverage_counter.get(False, {}).get(key, Counter())
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

def dump(self, read_data_provider, out_prefix):
    self.main_df = pd.DataFrame(self._calc_coverage())

    categories = { cat.id: cat for cat in read_data_provider.adm.get_categories() }

    gene_category_map = self._get_gene_category_map(self, categories, read_data_provider)

    coverage_columns = ["uniq_depth", "uniq_depth_covered", "uniq_horizontal", "combined_depth", "combined_depth_covered", "combined_horizontal"]

    for category in categories.values():
        features = pd.DataFrame.from_records(
            {"fid": feat.id, "feature": feat.name}
            for feat in self.adm.get_features(category=category.id)
        )

        cat_grouped = self._annotate_category_counts(
            self.main_df, 
            gene_category_map,
            coverage_columns,
            category.name
        ).groupby(category.name, as_index=False)
        
        coverage_df = cat_grouped[[category.name, "uniq_horizontal", "combined_horizontal",]].mean(numeric_only=True)
        depth_df = cat_grouped[[category.name, "uniq_depth", "uniq_depth_covered", "combined_depth", "combined_depth_covered",]].sum(numeric_only=True)
            
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
        
        new_order = ["feature", "uniq_depth", "uniq_depth_covered", "uniq_horizontal", "combined_depth", "combined_depth_covered", "combined_horizontal",]
        out_df[new_order].to_csv(
            f"{out_prefix}.{category.name}.coverage.txt", sep="\t", index=False, float_format="%.5f"
        )

    self.main_df.to_csv(out_prefix + ".all.coverage.txt", index=False, sep="\t", na_rep="NA")
    gene_category_map.to_csv(out_prefix + ".all.coverage_annotation.txt", index=False, sep="\t", na_rep="NA")

    read_data_provider.adm.dump(out_prefix + ".db")
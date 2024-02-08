import pandas as pd

class PandaProfiler:
	def __init__(self):
		self.main_df = None

	def dump(self, out_prefix):
		self.main_df.to_csv(
			f"{out_prefix}.panda_main_df.tsv",
			sep="\t",
			index=False,
			float_format="%.5f"
		)
	
	def add_records(self, hits):

		hits_df = pd.DataFrame(hits)
		hits_df["contrib"] = 1 / hits_df["n_aln"] / hits_df["library_mod"]
		hits_df["length"] = hits_df["end"] - hits_df["start"] + 1

		keep_columns = ["rid", "start", "end", "length", "contrib"]
		contrib_sums_uniq = hits_df[hits_df["is_ambiguous"] == False][keep_columns] \
			.groupby(by=keep_columns[:-1], as_index=False) \
			.sum(numeric_only=True)
		contrib_sums_combined = hits_df[keep_columns] \
			.groupby(by=keep_columns[:-1], as_index=False) \
			.sum(numeric_only=True)
		
		raw_counts_df = pd.merge(
			contrib_sums_uniq.drop(["length",], axis=1),
			contrib_sums_combined,
			on=("rid", "start", "end",),
			left_index=False, right_index=False,
			how="outer"
		) \
			.rename({"contrib_x": "uniq_raw", "contrib_y": "combined_raw"}, axis=1) \
			.fillna(0)
		
		if self.main_df is None:
			self.main_df = raw_counts_df
		else:
			self.main_df = pd.concat(
				(self.main_df, raw_counts_df,)
			).groupby(by=["gene", "length"], as_index=False).sum(numeric_only=True)
	

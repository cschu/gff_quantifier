import pandas as pd

class PandaProfiler:
	def __init__(self):
		self.main_df = None

	def get_gene_coords(self):
		for rid, start, end in zip(
			self.main_df["rid"], self.main_df["start"], self.main_df["end"]
		):
			yield rid, start, end

	def _annotate_records(self, gene_coords, refmgr, seqdb):
		gene_df = pd.DataFrame.from_records(
            { 
                "rid": rid,
                "start": start,
                "end": end,
                "gene": seqdb.get_db_sequence(
                    refmgr.get(rid[0] if isinstance(rid, tuple) else rid)[0],
                    start=start, end=end
                )[0].featureid,
			}
            # for rid, start, end in zip(df["rid"], df["start"], df["end"])
			for rid, start, end in gene_coords
		) \
			.drop_duplicates(keep="first")
		# gene_df.to_csv(self.out_prefix + ".gene_d.tsv", sep="\t", index=False)
		# hit_cols = ["gene", "rid", "start", "end", "rev_strand", "cov_start", "cov_end", "has_annotation", "n_aln", "is_ambiguous", "mate_id", "library_mod"]
		
		self.main_df = pd.merge(
			self.main_df,
			gene_df,
			on=("rid", "start", "end",),
			left_index=False, right_index=False,
			how="inner",
		) #[hit_cols]

	def profile(self, read_data_provider):
		self._annotate_records(
			self.get_gene_coords(),
			read_data_provider.reference_manager,
			read_data_provider.adm,
		)
		



	def dump(self, out_prefix):
		self.main_df.to_csv(
			f"{out_prefix}.panda_main_df.tsv",
			sep="\t",
			index=False,
			float_format="%.5f"
		)
	
	def add_records(self, hits):

		# [2024-02-08 14:51:17,846] count_stream: 
		# (
		# ([4308 1       447     True    7       157     True    None    True    2       1], 5),
		# ([19050 1       834     True    1       148     True    None    True    2       2], 5),
		# ([13361 1       501     True    61      211     True    None    True    2       1], 5),
		# ([37306 3       581     False   3       115     True    None    True    2       1], 5),
		# ([18264 2       331     True    251     331     True    None    True    2       1], 5),
		# ([19050 1       834     False   1       80      True    None    True    2       1], 5),
		# ([13361 1       501     False   1       143     True    None    True    2       2], 5),
		# ([18264 2       331     False   183     331     True    None    True    2       2], 5),
		# ([37306 3       581     True    33      183     True    None    True    2       2], 5),
		# ([4308  1       447     False   1       89      True    None    True    2       2], 5))

		hits_df = pd.DataFrame(hits)
		hits_df["contrib"] = 1 / hits_df["n_aln"] / hits_df["library_mod"]
		# hits_df["length"] = hits_df["end"] - hits_df["start"] + 1

		keep_columns = ["rid", "start", "end", "contrib"]
		contrib_sums_uniq = hits_df[hits_df["is_ambiguous"] == False][keep_columns] \
			.groupby(by=keep_columns[:-1], as_index=False) \
			.sum(numeric_only=True)
		contrib_sums_combined = hits_df[keep_columns] \
			.groupby(by=keep_columns[:-1], as_index=False) \
			.sum(numeric_only=True)
		
		raw_counts_df = pd.merge(
			contrib_sums_uniq,
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
			) \
				.groupby(by=["rid", "start", "end"], as_index=False) \
				.sum(numeric_only=True)
	

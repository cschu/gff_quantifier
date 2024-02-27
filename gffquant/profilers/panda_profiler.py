import pandas as pd

class PandaProfiler:
	def __init__(self):
		self.main_df = None

	def get_gene_coords(self):
		for rid, start, end in zip(
			self.main_df["rid"], self.main_df["start"], self.main_df["end"]
		):
			yield rid, start, end

	def _get_gene_category_map(self, categories, read_data_provider):
		gene_category_map = pd.merge(
			pd.DataFrame.from_records(
				self._get_gene_annotation(
					self.main_df[["rid", "start", "end"]],
					categories,
					read_data_provider.reference_manager,
					read_data_provider.adm,
				)
			),
			self.main_df[["gene", "rid", "start", "end"]],
			left_on=("refid", "start", "end",),
			right_on=("rid", "start", "end",),
			left_index=False, right_index=False,
			how="inner",
		)

		gene_category_map.to_csv(
			f"{read_data_provider.out_prefix}.gcmap.tsv",
			sep="\t",
			index=False,
		)
		return gene_category_map

	def _annotate_category_counts(self, counts_df, annotation_df, columns, category) -> pd.DataFrame:
		return pd.merge(
            annotation_df[["refid", "start", "end", "refname", category]],
            counts_df,
            left_index=False, right_index=False,
            left_on=("refname",),
            right_on=("gene",),
        ) \
        .dropna(axis=0, subset=[category,], how="any") #\
        # .explode(category, ignore_index=True)[[category,] + columns]

	def _get_gene_annotation(self, df, categories, refmgr, dbseq):
		for rid, start, end in zip(df["rid"], df["start"], df["end"]):
			ref, _ = refmgr.get(rid)
			for annseq in dbseq.get_db_sequence(ref, start=start, end=end):
				if annseq.annotation_str is not None: #and start == annseq.start and annseq.end == end:
					d = {"refid": rid, "start": start, "end": end, "refname": annseq.featureid}
					d.update({cat.name: None for cat in categories.values()})
					# annotated_cols.append(d)
					for item in annseq.annotation_str.split(";"):
						catid, features = item.split("=")
						d[categories.get(int(catid)).name] = [int(feat) for feat in features.split(",")]
					yield d

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

		# categories = { cat.id: cat for cat in seqdb.get_categories() }
		
		self.main_df = pd.merge(
			self.main_df,
			gene_df,
			on=("rid", "start", "end",),
			left_index=False, right_index=False,
			how="inner",
		) #[hit_cols]
		self.main_df["length"] = (self.main_df["end"] - self.main_df["start"] + 1)


	def profile(self, read_data_provider):
		self._annotate_records(
			self.get_gene_coords(),
			read_data_provider.reference_manager,
			read_data_provider.adm,
		)
		
		self.main_df["uniq_lnorm"] = self.main_df["uniq_raw"] / self.main_df["length"]
		self.main_df["combined_lnorm"] = self.main_df["combined_raw"] / self.main_df["length"]
		categories = { cat.id: cat for cat in read_data_provider.adm.get_categories() }

		gene_category_map = self._get_gene_category_map(categories, read_data_provider)

		count_columns = ["uniq_raw", "combined_raw", "uniq_lnorm", "combined_lnorm"]
		for category in categories.values():
			features = pd.DataFrame.from_records(
				{"fid": feat.id, "feature": feat.name }
				for feat in read_data_provider.adm.get_features(category=category.id)
			)
			category_counts = self._annotate_category_counts(
				self.main_df,
				gene_category_map,
				count_columns,
				category.name,
			) 
			
			category_counts.to_csv(
				f"{read_data_provider.out_prefix}.unexploded.{category.name}.pd.txt",
				sep="\t",
				index=False,
				float_format="%.5f",
			)

			category_row = [
				category_counts[col].sum(numeric_only=True)
				for col in ["uniq_raw", "uniq_lnorm", "combined_raw", "combined_lnorm",]
			]
			
			category_counts = category_counts \
				.explode(category.name, ignore_index=True)[[category.name,] + count_columns] \
				.groupby(category.name, as_index=False) \
				.sum(numeric_only=True)
			
			cat_df = pd.merge(
				features,
				category_counts,
				left_index=False, right_index=False,
				left_on=("fid",),
				right_on=(category.name,),
			) \
				.drop([category.name, "fid",], axis=1) \
				.sort_values(by=["feature",])
			
			cat_df["uniq_scaled"] = cat_df["uniq_lnorm"] * (
				cat_df["uniq_raw"].sum(numeric_only=True) / cat_df["uniq_lnorm"].sum(numeric_only=True)
			)
			cat_df["combined_scaled"] = cat_df["combined_lnorm"] * (
				cat_df["combined_raw"].sum(numeric_only=True) / cat_df["combined_lnorm"].sum(numeric_only=True)
			)
			out_cols = ["feature", "uniq_raw", "uniq_lnorm", "uniq_scaled", "combined_raw", "combined_lnorm", "combined_scaled"]

			header_rows = pd.DataFrame(
				[
					["total_reads"] + [read_data_provider.aln_counter["read_count"]] * (len(out_cols) - 1),
					["filtered_reads"] + [read_data_provider.aln_counter["filtered_read_count"]] * (len(out_cols) - 1),
					["category"] + category_row[:2] + [category_row[1] * (category_row[0] / category_row[1]),] + category_row[2:] + [category_row[3] * (category_row[2] / category_row[3]),],
				],
				columns = out_cols
			)
			


			pd.concat(
				[
					header_rows,
					cat_df[out_cols],
				]
			) \
				.to_csv(
					f"{read_data_provider.out_prefix}.{category.name}.pd.txt",
					sep="\t",
					index=False,
					float_format="%.5f",
				)


	def dump(self, out_prefix):
		
		self.main_df.to_csv(
			f"{out_prefix}.panda_main_df.tsv",
			sep="\t",
			index=False,
			float_format="%.5f"
		)

		out_cols = ["gene", "uniq_raw", "uniq_lnorm", "uniq_scaled", "combined_raw", "combined_lnorm", "combined_scaled"]

		self.main_df["uniq_scaled"] = self.main_df["uniq_lnorm"] * (
			self.main_df["uniq_raw"].sum(numeric_only=True) / self.main_df["uniq_lnorm"].sum(numeric_only=True)
		)
		self.main_df["combined_scaled"] = self.main_df["combined_lnorm"] * (
			self.main_df["combined_raw"].sum(numeric_only=True) / self.main_df["combined_lnorm"].sum(numeric_only=True)
		)

		self.main_df[out_cols] \
			.sort_values(by=["gene",]) \
			.to_csv(
				f"{out_prefix}.gene_counts.pd.txt",
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
	

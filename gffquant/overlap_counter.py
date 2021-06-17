import time
from collections import Counter

import numpy as np


"""
normalizeCounts nmethod counts sizes
	| nmethod `elem` [NMScaled, NMFpkm] = do
		-- count vectors always include a -1 at this point (it is
		-- ignored in output if the user does not request it, but is
		-- always computed). Thus, we compute the sum without it and do
		-- not normalize it later:
		let totalCounts v = withVector v (VU.sum . VU.tail)
		initial <- totalCounts counts
		normalizeCounts NMNormed counts sizes
		afternorm <- totalCounts counts
		let factor
				| nmethod == NMScaled = initial / afternorm
				| otherwise = 1.0e9 / initial --- 1e6 [million fragments] * 1e3 [kilo basepairs] = 1e9
		liftIO $ forM_ [1.. VUM.length counts - 1] (VUM.unsafeModify counts (* factor))
"""

class OverlapCounter(dict):
	COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled"]

	@staticmethod
	def normalise_counts(counts, feature_len, scaling_factor):
		'''Returns raw, length-normalised, and scaled feature counts.'''
		normalised = counts / feature_len
		scaled = normalised * scaling_factor
		return counts, normalised, scaled
	def __init__(self, out_prefix, db, do_overlap_detection=True, strand_specific=False):
		self.out_prefix = out_prefix
		self.db = db
		self.seqcounts = Counter()
		self.ambig_seqcounts = Counter()
		self.featcounts = dict()
		self.unannotated_reads = 0
		self.ambig_counts = dict()
		self.has_ambig_counts = False
		self.feature_scaling_factors = dict()
		self.do_overlap_detection = do_overlap_detection
		self.strand_specific = strand_specific
		self.gene_counts = dict()
		#self.ambig_gene_counts = Counter()

	def update_unique_counts(self, rid, aln=None, rev_strand=False):
		'''
		Generates a hit list from the overlaps resulting from an intervaltree query,
		adds the number of alternative alignments and stores the results for each reference sequence.
		'''
		if aln:
			overlaps = self.db.get_overlaps(*aln)
			if overlaps:
				self.setdefault(rid, Counter()).update((ovl.begin, ovl.end, rev_strand) for ovl in overlaps)
			else:
				self.unannotated_reads += 1
		if self.strand_specific and not self.do_overlap_detection:
			self.seqcounts[(rid, rev_strand)] += 1
		else:
			self.seqcounts[rid] += 1

	def _compute_count_vector(self, bins, rid, aln, strand, strand_specific):
		# calculate the count vector for the current region
		# count vectors are of the format (raw_uniq, normed_uniq, raw_ambi, normed_ambi)
		# att: _ambi counts are unique + ambiguous
		counts = np.zeros(bins)
		region_length = aln[1] - aln[0] + 1
		counts[0] = counts[1] = self.get(rid, dict()).get(aln, 0.0)
		counts[1] /= region_length
		counts[2] = counts[3] = counts[0] + self.ambig_counts.get(rid, dict()).get(aln, 0.0)
		counts[3] /= region_length
		if strand_specific:
			rev_strand = aln[-1]
			is_antisense = (strand == "+" and rev_strand) or (strand == "-" and not rev_strand)
			bin_offset = 8 if is_antisense else 4
			counts[bin_offset:bin_offset + 4] += counts[:4]

		return counts

	def _compute_genes_count_vector(self, rid, length, strand_specific=False):
		counts = np.zeros(12 if strand_specific else 4)
		if strand_specific:
			PLUS_STRAND, MINUS_STRAND = True, False
			uniq_plus = self.seqcounts[(rid, PLUS_STRAND)]
			uniq_minus = self.seqcounts[(rid, MINUS_STRAND)]
			ambig_plus = self.ambig_seqcounts[(rid, PLUS_STRAND)]
			ambig_minus = self.ambig_seqcounts[(rid, MINUS_STRAND)]
			counts[0] = counts[1] = uniq_plus + uniq_minus
			counts[1] /= length
			counts[2] = counts[3] = counts[0] + ambig_plus + ambig_minus
			counts[3] /= length
			counts[4] = counts[5] = uniq_plus
			counts[5] /= length
			counts[6] = counts[7] = counts[4] + ambig_plus
			counts[7] /= length
			counts[8] = counts[9] = uniq_minus
			counts[9] /= length
			counts[10] = counts[11] = counts[8] + ambig_minus
			counts[11] /= length
		else:
			counts[0] = counts[1] = self.seqcounts[rid]
			counts[1] /= length
			counts[2] = counts[3] = counts[0] + self.ambig_seqcounts[rid]
			counts[3] /= length
		return counts

	def _iterate_database(self, bins, bam, strand_specific):
		total_counts, feature_count_sums = np.zeros(4), dict()

		for ref, region_annotation in self.db.iterate():
			rid = bam.revlookup_reference(ref)
			if rid is not None:
				_, region_length = bam.get_reference(rid)
				counts = self._compute_genes_count_vector(rid, region_length, strand_specific=strand_specific)
				total_counts, feature_count_sums = self._distribute_feature_counts(
					bins, counts, region_annotation[1:], total_counts, feature_count_sums
				)
				self._add_count_vector(counts, ref, self.gene_counts, bins)

		return total_counts, feature_count_sums


	def _iterate_bedcounts(self, bins, feature_lengths, strand_specific):
		total_counts, feature_count_sums = np.zeros(4), dict()

		for rid in set(self.seqcounts).union(self.ambig_seqcounts):
			region_length = feature_lengths.get(rid)
			feature = rid.split("::")[-1]

			if region_length is None:
				raise ValueError(f"Cannot determine length of reference {rid}.")
			counts = self._compute_genes_count_vector(rid, region_length, strand_specific=strand_specific)

			total_counts, feature_count_sums = self._distribute_feature_counts(
				bins, counts, [("feature", (feature,))], total_counts, feature_count_sums
			)

		return total_counts, feature_count_sums


	def _distribute_feature_counts(
		self, bins, counts, region_annotation, total_counts, feature_count_sums
	):
		for ftype, ftype_counts in region_annotation:
			total_fcounts = feature_count_sums.setdefault(ftype, np.zeros(4))
			for i, ft_ct in enumerate(ftype_counts, start=1):
				fcounts = self.featcounts.setdefault(ftype, dict())
				self._add_count_vector(counts, ft_ct, fcounts, bins)

			try:
				inc = counts[:4] * i
			except NameError:
				inc = 0
			total_counts += inc
			total_fcounts += inc

		return total_counts, feature_count_sums

	def _add_count_vector(self, count_vector, target_id, target_ctr, bins):
		counts = target_ctr.setdefault(target_id, np.zeros(bins))
		counts += count_vector

	def _iterate_counts(self, bins, bam, strand_specific):
		total_counts, feature_count_sums = np.zeros(4), dict()

		for rid in set(self.keys()).union(self.ambig_counts):
			ref = bam.get_reference(rid)[0]
			regions = set(self.get(rid, set())).union(self.ambig_counts.get(rid, set()))
			for start, end, rev_strand in regions:
				# region_annotation is a tuple of key-value pairs:
				# (strand, func_category1: subcategories, func_category2: subcategories, ...)
				# the first is the strand, the second is the gene id, the rest are the features
				region_annotation = self.db.get_data(ref, start, end)
				counts = self._compute_count_vector(
					bins, rid, (start, end, rev_strand), region_annotation[0][1], strand_specific
				)
				# how to extract the gene counts in genome mode?
				self._add_count_vector(counts, region_annotation[1][1], self.gene_counts, bins)

				# distribute the counts to the associated functional (sub-)categories
				total_counts, feature_count_sums = self._distribute_feature_counts(
					bins, counts, region_annotation[2:], total_counts, feature_count_sums
				)

		return total_counts, feature_count_sums


	def annotate_counts(self, bamfile=None, feature_lengths=None, itermode="counts"):
		"""
		Distributes read counts against aligned regions to the associated functional categories.
		"""
		print("Processing counts ...", flush=True)
		t0 = time.time()
		bins = 12 if self.strand_specific else 4


		annotation_f = {
			"counts": self._iterate_counts,
			"database": self._iterate_database,
			"bedcounts": self._iterate_bedcounts
		}.get(itermode)

		if not annotation_f:
			raise ValueError(f"Unknown annotation_mode {itermode}")

		if bool(bamfile) == bool(feature_lengths):
			raise ValueError("Cannot determine feature lengths.")
		if not feature_lengths:
			feature_lengths = bamfile

		total_counts, feature_count_sums = annotation_f(bins, feature_lengths, self.strand_specific)

		# calculate the scaling factors
		total, total_normed, total_ambi, total_ambi_normed = total_counts

		default_scaling_factor = 0
		self.feature_scaling_factors = {
			"total": (total / total_normed) if total_normed else default_scaling_factor,
			"total_ambi": (total_ambi / total_ambi_normed) if total_ambi_normed else default_scaling_factor
		}

		for ftype, counts in feature_count_sums.items():
			total, total_normed, total_ambi, total_ambi_normed = counts
			self.feature_scaling_factors[ftype] = (
				(total / total_normed) if total_normed else default_scaling_factor,
				(total_ambi / total_ambi_normed) if total_ambi_normed else default_scaling_factor
			)

		self.clear()
		self.ambig_counts.clear()
		t1 = time.time()
		print("Processed counts in {n_seconds}s.".format(n_seconds=t1-t0), flush=True)


	def update_ambiguous_counts(self, hits, n_aln, unannotated=0, feat_distmode="all1"):
		self.has_ambig_counts = True
		self.unannotated_reads += unannotated
		if self.strand_specific and not self.do_overlap_detection:
			n_total = sum(self.seqcounts[(rid, True)] + self.seqcounts[(rid, False)] for rid in hits)
		else:
			n_total = sum(self.seqcounts[rid] for rid in hits)

		increment = (1 / n_aln) if feat_distmode == "1overN" else 1
		for rid, regions in hits.items():
			if self.do_overlap_detection:
				for start, end, rev_str in regions:
					self.ambig_counts.setdefault(rid, Counter())[(start, end, rev_str)] += increment
#				reg_count = self.ambig_counts.setdefault(rid, Counter())
#
#				if feat_distmode == "all1":
#					increment = 1
#				elif feat_distmode == "1overN":
#					increment = 1 / n_aln
#				else:
#					uniq_counts = self.get((start, end, True), 0) + self.get((start, end, False), 0)
#					if n_total and uniq_counts:
#						increment = uniq_counts /
#
#				reg_count[(start, end, rev_str)] += increment

#290			 hits.setdefault(rid, set()).add((start, end, SamFlags.is_reverse_strand(flag)))

				if n_total and self.seqcounts[rid]:
					self.ambig_seqcounts[rid] += self.seqcounts[rid] / n_total * len(hits)
				else:
					self.ambig_seqcounts[rid] += 1 / len(hits)
			else:
				start, end, rev_str = list(regions)[0]
				key = (rid, rev_str) if self.strand_specific else rid
				self.ambig_seqcounts[key] += (1 / n_aln) if feat_distmode == "1overN" else 1

	@staticmethod
	def calculate_seqcount_scaling_factor(counts, bam):
		raw_total = sum(counts.values())
		normed_total = sum(count / bam.get_reference(rid)[1] for rid, count in counts.items())
		return raw_total / normed_total
	@staticmethod
	def calculate_feature_scaling_factor(counts, include_ambig=False):
		raw_total, normed_total = 0, 0
		for raw, norm, raw_ambi, norm_ambi, *_ in counts.values():
			raw_total += raw
			normed_total += norm
			if include_ambig:
				raw_total += raw_ambi
				normed_total += norm_ambi
		return raw_total / normed_total


	def _compile_output_row(self, counts, scaling_factor=1, ambig_scaling_factor=1):
		p, row = 0, list()
		# unique counts
		row.extend(counts[p:p + 2])
		row.append(row[-1] * scaling_factor)
		p += 2
		# ambiguous counts
		if self.has_ambig_counts:
			row.extend(counts[p:p + 2])
			row.append(row[-1] * ambig_scaling_factor)
			p += 2
		# sense-strand unique
		if self.strand_specific:
			row.extend(counts[p:p + 2])
			row.append(row[-1] * scaling_factor)
			p += 2
			# sense-strand ambiguous
			if self.has_ambig_counts:
				row.extend(counts[p:p + 2])
				row.append(row[-1] * ambig_scaling_factor)
				p += 2
			# antisense-strand unique
			row.extend(counts[p:p + 2])
			row.append(row[-1] * scaling_factor)
			p += 2
			if self.has_ambig_counts:
				row.extend(counts[p:p + 2])
				row.append(row[-1] * ambig_scaling_factor)
		return row

	def get_header(self):
		header = list()
		header.extend(
			f"uniq_{element}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS
		)
		if self.has_ambig_counts:
			header.extend(
				f"combined_{element}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS
			)
		if self.strand_specific:
			for strand in ("ss", "as"):
				header.extend(
					f"uniq_{element}_{strand}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS
				)
				if self.has_ambig_counts:
					header.extend(
						f"combined_{element}_{strand}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS
					)
		return header


	def _dump_feature_counts(self):
		with open(f"{self.out_prefix}.feature_counts.txt", "w") as feat_out:
			print("subfeature", *self.get_header(), sep="\t", file=feat_out, flush=True)
			print("unannotated", self.unannotated_reads, sep="\t", file=feat_out, flush=True)
			for ftype, counts in sorted(self.featcounts.items()):
				print(f"#{ftype}", file=feat_out, flush=True)
				scaling_factor, ambig_scaling_factor = self.feature_scaling_factors[ftype]
				for subf, sf_counts in sorted(counts.items()):
					out_row = self._compile_output_row(
						sf_counts, scaling_factor=scaling_factor, ambig_scaling_factor=ambig_scaling_factor
					)
					print(
						subf, out_row[0], *(f"{c:.5f}" for c in out_row[1:]),
						flush=True, sep="\t", file=feat_out
					)


	def _dump_seq_counts(self, bam):
		SEQ_COUNT_HEADER = ["seqid_int", "seqid", "length"] + OverlapCounter.COUNT_HEADER_ELEMENTS

		if self.strand_specific and not self.do_overlap_detection:
			_seqcounts = Counter()
			for (rid, rev_strand), count in self.seqcounts.items():
				_seqcounts[rid] += count
			self.seqcounts = _seqcounts
			if self.has_ambig_counts:
				_seqcounts = Counter()
				for (rid, rev_strand), count in self.ambig_seqcounts.items():
					_seqcounts[rid] += count
				self.ambig_seqcounts = _seqcounts

		with open(f"{self.out_prefix}.seqname.uniq.txt", "w") as seq_out:
			print(*SEQ_COUNT_HEADER, sep="\t", flush=True, file=seq_out)
			if sum(self.seqcounts.values()):
				seqcount_scaling_factor = OverlapCounter.calculate_seqcount_scaling_factor(self.seqcounts, bam)
				for rid, count in self.seqcounts.items():
					seq_id, seq_len = bam.get_reference(rid)
					norm, scaled = OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)[1:]
					print(
						rid, seq_id, seq_len, count, f"{norm:.5f}", f"{scaled:.5f}",
						flush=True, sep="\t", file=seq_out
					)

		if self.ambig_seqcounts:
			with open(f"{self.out_prefix}.seqname.dist1.txt", "w") as seq_out:
				print(*SEQ_COUNT_HEADER, sep="\t", flush=True, file=seq_out)
				self.seqcounts.update(self.ambig_seqcounts)
				seqcount_scaling_factor = OverlapCounter.calculate_seqcount_scaling_factor(self.seqcounts, bam)
				for rid, count in self.seqcounts.items():
					seq_id, seq_len = bam.get_reference(rid)
					raw, norm, scaled = OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)
					print(
						rid, seq_id, seq_len, f"{raw:.5f}", f"{norm:.5f}", f"{scaled:.5f}",
						flush=True, sep="\t", file=seq_out
					)


	def _dump_gene_counts(self, bam=None):
		with open(f"{self.out_prefix}.gene_counts.txt", "w") as gene_out:
			print("gene", *self.get_header(), sep="\t", file=gene_out, flush=True)
			if self.do_overlap_detection:
				gene_counts = self.gene_counts
				scaling_factor = self.feature_scaling_factors["total"]
				ambig_scaling_factor = self.feature_scaling_factors["total_ambi"]
			else:
				gene_counts = dict()
				gene_ids = set(
					(gene_id[0] if isinstance(gene_id, tuple) else gene_id)
					for gene_id in set(self.seqcounts).union(self.ambig_seqcounts)
				)
				gene_counts = {
					gene_id: self._compute_genes_count_vector(
						gene_id, bam.get_reference(gene_id)[1], strand_specific=self.strand_specific
					)
					for gene_id in gene_ids
				}
				counts_for_scaling = np.zeros(4)
				for counts in gene_counts.values():
					counts_for_scaling += counts[:4]
				scaling_factor = counts_for_scaling[0] / counts_for_scaling[1]
				ambig_scaling_factor = counts_for_scaling[2] / counts_for_scaling[3]

			for gene, g_counts in sorted(gene_counts.items()):
				out_row = self._compile_output_row(
					g_counts,
					scaling_factor=scaling_factor,
					ambig_scaling_factor=ambig_scaling_factor
				)
				if isinstance(gene, tuple):
					gene = gene[0]
				if not self.do_overlap_detection:
					gene, _ = bam.get_reference(gene)
				print(gene, out_row[0], *(f"{c:.5f}" for c in out_row[1:]), flush=True, sep="\t", file=gene_out)


	def dump_counts(self, bam=None):
		print("Dumping overlap counters...", flush=True)
		print("Has ambiguous counts:", self.has_ambig_counts, flush=True)

		self._dump_feature_counts()

		if self.gene_counts:
			self._dump_gene_counts(bam=bam)

		if bam:
			self._dump_seq_counts(bam)

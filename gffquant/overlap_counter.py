import os
import time
import sys
from collections import Counter

import numpy


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
	@staticmethod
	def normalise_counts(counts, feature_len, scaling_factor):
		'''Returns raw, length-normalised, and scaled feature counts.'''
		normalised = counts / feature_len
		scaled = normalised * scaling_factor
		return counts, normalised, scaled
	def __init__(self, out_prefix):
		self.out_prefix = out_prefix
		self.seqcounts = Counter()
		self.ambig_seqcounts = Counter()
		self.featcounts = dict()
		self.unannotated_reads = 0
		self.ambig_counts = dict()
		self.has_ambig_counts = False
		self.scaling_factor = None
		self.scaling_factor_ambi = None
		pass
	def update_unique_counts(self, rid, overlaps, rev_strand=False):
		'''
		Generates a hit list from the overlaps resulting from an intervaltree query,
		adds the number of alternative alignments and stores the results for each reference sequence.
		'''
		if overlaps:
			self.setdefault(rid, Counter()).update((ovl.begin, ovl.end, rev_strand) for ovl in overlaps)
		else:
			self.unannotated_reads += 1
		self.seqcounts[rid] += 1

	def annotate_counts(self, bam, gff_dbm, strand_specific=False):
		print("Processing counts ...", flush=True)
		t0 = time.time()

		total, total_ambi = 0, 0
		total_normed, total_ambi_normed = 0, 0

		for rid in set(self.keys()).union(self.ambig_counts):
			ref = bam.get_reference(rid)[0]
			for start, end, rev_strand in set(self.get(rid, set())).union(self.ambig_counts.get(rid, set())):
				region_length = end - start + 1
				region_annotation = gff_dbm.get_data(ref, start, end)
				strand = region_annotation[0][1]

				is_antisense = strand_specific and ((strand == "+" and rev_strand) or (strand == "-" and not rev_strand))
				if strand_specific:
					bins, bin_offset = 12, (8 if is_antisense else 4)
				else:
					bins, bin_offset = 4, 0

				counts = [self.get(rid, dict()).get((start, end, rev_strand), 0), 0, 0, 0]
				counts[1] = counts[0] / region_length
				counts[2] = counts[0] + self.ambig_counts.get(rid, dict()).get((start, end, rev_strand), 0)
				counts[3] = counts[2] / region_length

				total += counts[0]
				total_normed += counts[1]
				total_ambi += counts[2]
				total_ambi_normed += counts[3]

				for ftype, ftype_counts in region_annotation[1:]:
					for ft_ct in ftype_counts:
						fcounts = self.featcounts.setdefault(ftype, dict()).setdefault(ft_ct, numpy.array([0.0 for i in range(bins)]))
						for i, c in enumerate(counts):
							fcounts[i] += c
							if bin_offset:
								fcounts[i + bin_offset] += c

				#for ftype, values in region_annotation[1:]:
				#	for v in values:
				#		#count = self.get(rid, dict()).get((start, end, rev_strand), 0)
				#		#ambi_count = count + self.ambig_counts.get(rid, dict()).get((start, end, rev_strand), 0)
				#		#if strand_specific:
				#		#	bins, bin_offset = 12, (8 if is_antisense else 4)
				#		#else:
				#		#	bins, bin_offset = 4, 0
				#		fcounts = self.featcounts.setdefault(ftype, dict()).setdefault(v, numpy.array([0.0 for i in range(bins)]))
				#		for i, c in enumerate((count, count / region_length, ambi_count, ambi_count / region_length)):
				#			fcounts[i] += c
				#			if bin_offset:
				#				fcounts[i + bin_offset] += c

		self.scaling_factor = (total / total_normed) if total_normed else None
		self.scaling_factor_ambi = (total_ambi / total_ambi_normed) if total_ambi_normed else None

		self.clear()
		self.ambig_counts.clear()
		t1 = time.time()
		print("Processed counts in {n_seconds}s.".format(n_seconds=t1-t0), flush=True)







	def update_ambiguous_counts(self, hits, n_aln, unannotated, gff_dbm, bam, feat_distmode="all1", strand_specific=False):
		self.has_ambig_counts = True
		n_total = sum(self.seqcounts[rid] for rid in hits)
		for rid, regions in hits.items():
			for start, end, rev_str in regions:
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
				
				self.ambig_counts.setdefault(rid, Counter())[(start, end, rev_str)] += (1 / n_aln) if feat_distmode == "1overN" else 1

			if n_total and self.seqcounts[rid]:
				self.ambig_seqcounts[rid] += self.seqcounts[rid] / n_total * len(hits)
			else:
				self.ambig_seqcounts[rid] += 1 / len(hits)

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

	def dump_counts(self, bam, strand_specific=False):

		COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled"]
		SEQ_COUNT_HEADER = ["seqid_int", "seqid", "length"] + COUNT_HEADER_ELEMENTS
		counts_template = "{:.5f}\t{:.5f}\t{:.5f}"

		print("Dumping overlap counters...", flush=True)
		print("Has ambiguous counts:", self.has_ambig_counts, flush=True)
		with open("{prefix}.seqname.uniq.txt".format(prefix=self.out_prefix), "w") as seq_out:
			print(*SEQ_COUNT_HEADER, sep="\t", flush=True, file=seq_out)
			seqcount_scaling_factor = OverlapCounter.calculate_seqcount_scaling_factor(self.seqcounts, bam)
			for rid, count in self.seqcounts.items():
				seq_id, seq_len = bam.get_reference(rid)
				print(rid, seq_id, seq_len, count, "{:.5f}\t{:.5f}".format(*OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)[1:]), flush=True, sep="\t", file=seq_out)

		with open("{prefix}.feature_counts.txt".format(prefix=self.out_prefix), "w") as feat_out:

			header = ["subfeature"]
			header.extend("uniq_{}".format(element) for element in COUNT_HEADER_ELEMENTS)
			if self.has_ambig_counts:
				header.extend("ambig_{}".format(element) for element in COUNT_HEADER_ELEMENTS)
			if strand_specific:
				for strand in ("ss", "as"):
					header.extend("uniq_{}_{}".format(element, strand) for element in COUNT_HEADER_ELEMENTS)
					if self.has_ambig_counts:
						header.extend("ambig_{}_{}".format(element, strand) for element in COUNT_HEADER_ELEMENTS)
			print(*header, sep="\t", file=feat_out, flush=True)

			print("unannotated", self.unannotated_reads, sep="\t", file=feat_out, flush=True)
			for ftype, counts in sorted(self.featcounts.items()):
				print("#{}".format(ftype), file=feat_out, flush=True)
				#feature_scaling_factor = OverlapCounter.calculate_feature_scaling_factor(counts)
				#feature_scaling_factor_ambig = OverlapCounter.calculate_feature_scaling_factor(counts, include_ambig=True)

				for subf, sf_counts in sorted(counts.items()):
					# first batch: unique
					out_row = list(sf_counts[:2])
					out_row.append(out_row[-1] * self.scaling_factor)
					# next batch: ambiguous (if exist)
					if self.has_ambig_counts:
						out_row.extend(sf_counts[2:4])
						out_row.append(out_row[-1] * self.scaling_factor_ambi)
					# next batch: sense-strand unique
					if strand_specific:
						out_row.extend(sf_counts[4:6])
						out_row.append(out_row[-1] * self.scaling_factor)
						# next batch: sense-strand ambiguous
						if self.has_ambig_counts:
							out_row.extend(sf_counts[6:8])
							out_row.append(out_row[-1] * self.scaling_factor_ambi)
						# next batch antisense-strand unique
						out_row.extend(sf_counts[8:10])
						out_row.append(out_row[-1] * self.scaling_factor)
						# next batch: antisense-strand ambiguous
						if self.has_ambig_counts:
							out_row.extend(sf_counts[10:12])
							out_row.append(out_row[-1] * self.scaling_factor_ambi)

					print(subf, out_row[0], *("{:.5f}".format(c) for c in out_row[1:]), flush=True, sep="\t", file=feat_out)

		if self.ambig_seqcounts:
			with open("{prefix}.seqname.dist1.txt".format(prefix=self.out_prefix), "w") as seq_out:
				print(*SEQ_COUNT_HEADER, sep="\t", flush=True, file=seq_out)
				self.seqcounts.update(self.ambig_seqcounts)
				seqcount_scaling_factor = OverlapCounter.calculate_seqcount_scaling_factor(self.seqcounts, bam)
				for rid, count in self.seqcounts.items():
					seq_id, seq_len = bam.get_reference(rid)
					print(rid, seq_id, seq_len, counts_template.format(*OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)), flush=True, sep="\t", file=seq_out)

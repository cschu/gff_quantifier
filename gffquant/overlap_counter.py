import os
import time
import sys
from collections import Counter

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
		pass
	def update_unique_counts(self, rid, overlaps):
		'''
		Generates a hit list from the overlaps resulting from an intervaltree query,
		adds the number of alternative alignments and stores the results for each reference sequence.
		'''
		if overlaps:
			self.setdefault(rid, Counter()).update((ovl.begin, ovl.end) for ovl in overlaps)
		else:
			self.unannotated_reads += 1
		self.seqcounts[rid] += 1

	def annotate_counts(self, bam, gff_dbm):
		print("Processing counts ...", flush=True)
		t0 = time.time()
		for rid in set(self.keys()).union(self.ambig_counts):
			ref = bam.get_reference(rid)[0]
			for region in set(self.get(rid, set())).union(self.ambig_counts.get(rid, set())):
				region_length = region[1] - region[0] + 1
				region_annotation = gff_dbm.get_data(ref, *region)
				for ftype, values in region_annotation:
					for v in values:
						fcounts = self.featcounts.setdefault(ftype, dict()).setdefault(v, [0, 0, 0, 0])
						count = self.get(rid, dict()).get(region, 0)
						ambi_count = count + self.ambig_counts.get(rid, dict()).get(region, 0)
						fcounts[0] += count
						fcounts[1] += count / region_length
						fcounts[2] += ambi_count
						fcounts[3] += ambi_count / region_length
		self.clear()
		self.ambig_counts.clear()
		t1 = time.time()
		print("Processed counts in {n_seconds}s.".format(n_seconds=t1-t0), flush=True)

	def update_ambiguous_counts(self, hits, n_aln, unannotated, gff_dbm, bam, feat_distmode="all1"):
		if feat_distmode in ("all1", "1overN"):
			n_total = sum(self.seqcounts[rid] for rid in hits)
			for rid, regions in hits.items():
				for start, end in regions:
					self.ambig_counts.setdefault(rid, Counter())[(start, end)] += (1 / n_aln) if feat_distmode == "1overN" else 1
			# 	print("RID", rid, file=sys.stderr)

				try:
					self.ambig_seqcounts[rid] += self.seqcounts[rid] / n_total * len(hits)
				except ZeroDivisionError:
					self.ambig_seqcounts[rid] += 1 / len(hits)

	@staticmethod
	def calculate_seqcount_scaling_factor(counts, bam):
		raw_total = sum(counts.values())
		normed_total = sum(count / bam.get_reference(rid)[1] for rid, count in counts.items())
		return raw_total / normed_total
	@staticmethod
	def calculate_feature_scaling_factor(counts, include_ambig=False):
		raw_total, normed_total = 0, 0
		for raw, norm, raw_ambi, norm_ambi in counts.values():
			if include_ambig:
				raw_total += raw_ambi
				normed_total += norm_ambi
			else:
				raw_total += raw
				normed_total += norm
		return raw_total / normed_total

	def dump_counts(self, bam):
		counts_template = "{:.5f}\t{:.5f}\t{:.5f}"

		print("Dumping overlap counters...", flush=True)
		with open("{prefix}.seqname.uniq.txt".format(prefix=self.out_prefix), "w") as seq_out:
			seqcount_scaling_factor = OverlapCounter.calculate_seqcount_scaling_factor(self.seqcounts, bam)
			for rid, count in self.seqcounts.items():
				seq_id, seq_len = bam.get_reference(rid)
				print(rid, seq_id, seq_len, counts_template.format(*OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)), flush=True, sep="\t", file=seq_out)
		if self.ambig_seqcounts:
			with open("{prefix}.seqname.dist1.txt".format(prefix=self.out_prefix), "w") as seq_out:
				self.seqcounts.update(self.ambig_seqcounts)
				seqcount_scaling_factor = OverlapCounter.calculate_seqcount_scaling_factor(self.seqcounts, bam)
				for rid, count in self.seqcounts.items():
					seq_id, seq_len = bam.get_reference(rid)
					print(rid, seq_id, seq_len, counts_template.format(*OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)), flush=True, sep="\t", file=seq_out)

		with open("{prefix}.feature_counts.txt".format(prefix=self.out_prefix), "w") as feat_out:
			print("unannotated", self.unannotated_reads, sep="\t", file=feat_out, flush=True)
			for ftype, counts in sorted(self.featcounts.items()):
				print("#{}".format(ftype), file=feat_out, flush=True)
				feature_scaling_factor = OverlapCounter.calculate_feature_scaling_factor(counts)
				feature_scaling_factor_ambig = OverlapCounter.calculate_feature_scaling_factor(counts, include_ambig=True)

				for subf, (raw, norm, raw_ambi, norm_ambi) in sorted(counts.items()):
					print(
						subf,
						counts_template.format(raw, norm, norm * feature_scaling_factor),
						counts_template.format(raw_ambi, norm_ambi, norm_ambi * feature_scaling_factor_ambig),
						flush=True, sep="\t", file=feat_out
					)

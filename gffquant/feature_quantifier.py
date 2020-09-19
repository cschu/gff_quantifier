import os
import time
import sys
from collections import Counter
import json
import yaml
from datetime import datetime
from functools import lru_cache

from intervaltree import IntervalTree

from gffquant.bamreader import BamFile

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

class GffDatabaseManager:
	def _read_index(self, f):
		db_index = dict()
		for line in open(f, "rt"):
			line = line.strip().split("\t")
			db_index.setdefault(line[0], list()).append(list(map(int, line[1:3])))
		return db_index
	def __init__(self, db, db_index):
		self.db_index = self._read_index(db_index)
		self.db = open(db, "rt")
		self.loaded_data = None
		self.loaded_ref = None
		self.interval_tree = None
	@lru_cache(maxsize=4096)
	def _read_data(self, ref_id, include_payload=False):
		gff_annotation = dict()
		for offset, size in self.db_index.get(ref_id, list()):
			self.db.seek(offset)
			for line in self.db.read(size).strip("\n").split("\n"):
				if not line.startswith("#"):
					line = line.strip().split("\t")
					features = dict()
					if include_payload:
						features = tuple((item.split("=")[0], tuple(sorted(item.split("=")[1].split(",")))) for item in line[8].strip().split(";") if not item.startswith("ID"))
					key = (line[0], int(line[3]), int(line[4]) + 1)
					gff_annotation[key] = features
		if not gff_annotation and not include_payload:
			print("WARNING: contig {contig} does not have an annotation in the index.".format(contig=ref_id), file=sys.stderr, flush=True)
		return gff_annotation
	@lru_cache(maxsize=4096)
	def _get_tree(self, ref, cache_data=False):
		return IntervalTree.from_tuples(sorted([key[1:] for key in self._read_data(ref, include_payload=cache_data)]))
	def get_data(self, ref, start, end):
		return self._read_data(ref, include_payload=True).get((ref, start, end), dict())
	def get_overlaps(self, ref, start, end, cache_data=False):
		return self._get_tree(ref, cache_data=cache_data)[start:end]


class FeatureQuantifier:
	def _read_count_config(self, config):
		default = {"multiple": self.default_multiple_alignments, "normalization": self.default_normalization}
		self._count_config = yaml.load(open(config), Loader=yaml.SafeLoader) if config is not None else default

	def __init__(self, gff_db, gff_index, count_config=None, multiple_alignments="unique_only", normalization="scaled"):
		self.gff_dbm = GffDatabaseManager(gff_db, gff_index)
		self.default_multiple_alignments = multiple_alignments
		self.default_normalization = normalization
		self._read_count_config(count_config)
		self.umap_cache = dict()

	def process_unique_cache(self, rid, ref):
		for qname, uniq_aln in self.umap_cache.items():
			n_aln = len(uniq_aln)
			if n_aln > 2:
				print("WARNING: more than two primary alignments for {qname} ({n}). Ignoring.".format(qname=qname, n=n_aln), file=sys.stderr, flush=True)
				continue

			start, end = uniq_aln[0][1:] if n_aln == 1 else BamFile.calculate_fragment_borders(*uniq_aln[0][1:], *uniq_aln[1][1:])
			self.overlap_counter.update_unique_counts(rid, self.gff_dbm.get_overlaps(ref, start, end))

		self.umap_cache.clear()

	def process_ambiguous_alignments(self, bam):
		t0 = time.time()
		current_group = None

		n_align = 0
		for aln_count, aln in bam.get_alignments(allow_multiple=True, allow_unique=False, disallowed_flags=0x800):
			if current_group is None or current_group.qname != aln.qname:
				if current_group is not None:
					n_align += current_group.n_align()
					current_group.resolve(self.overlap_counter, self.gff_dbm, bam)
					print("{n_align} alignments processed.".format(n_align=n_align), file=sys.stderr, flush=True)
				current_group = AlignmentGroup(aln)
			else:
				current_group.add_alignment(aln)

		if current_group is not None:
			n_align += current_group.n_align()
			current_group.resolve(self.overlap_counter, self.gff_dbm, bam)
			print("{n_align} alignments processed.".format(n_align=n_align), file=sys.stderr, flush=True)

		t1 = time.time()
		print("Processed {n_align} secondary alignments in {n_seconds:.3f}s.".format(
			n_align=aln_count, n_seconds=t1-t0), flush=True)


	def process_unique_alignments(self, bam):
		t0 = time.time()
		current_ref, current_rid = None, None

		for aln_count, aln in bam.get_alignments(allow_multiple=False, disallowed_flags=0x800):
			ref = bam.get_reference(aln.rid)[0]
			start, end = aln.start, aln.end
			qname = aln.qname

			if ref != current_ref:
				print("{time}\tNew reference: {ref} ({rid}/{n_ref}). {n_aln} alignments processed.".format(
					time=datetime.now().strftime("%m/%d/%Y,%H:%M:%S"),
					ref=ref, rid=aln.rid, n_ref=bam.n_references(), n_aln=aln_count),
					file=sys.stderr, flush=True)

				self.process_unique_cache(current_rid, current_ref)
				current_rid, current_ref = aln.rid, ref

			if aln.is_paired() and aln.rid == aln.rnext and not aln.flag & 0x8:
				# need to keep track of read pairs to avoid dual counts
				# check if the mate has already been seen
				mates = self.umap_cache.setdefault(qname, list())
				if mates:
					# if it has, calculate the total fragment size and remove the pair from the cache
					if aln.rnext != mates[0][0]:
						print("WARNING: alignment {qname} seems to be corrupted: {aln1} {aln2}".format(qname=qname, aln1=str(aln), aln2=str(mates[0])), flush=True, file=sys.stderr)
						continue # i don't think this ever happens
					else:
						start, end = BamFile.calculate_fragment_borders(aln.start, aln.end, *mates[0][1:])
					del self.umap_cache[qname]
				else:
					# otherwise cache the first encountered mate and advance to the next read
					mates.append((aln.rid, aln.start, aln.end))
					continue

			# at this point only single-end reads and merged pairs should be processed here
			self.overlap_counter.update_unique_counts(aln.rid, self.gff_dbm.get_overlaps(current_ref, start, end))

		# clear cache
		self.process_unique_cache(current_rid, current_ref)
		t1 = time.time()

		print("Processed {n_align} primary alignments in {n_seconds:.3f}s.".format(
			n_align=aln_count, n_seconds=t1-t0), flush=True
		)

	def process_data(self, bamfile, bamfile_ns=None, out_prefix="gffquant"):
		self.overlap_counter = OverlapCounter(out_prefix)
		bam = BamFile(bamfile)

		self.process_unique_alignments(bam)

		if bamfile_ns is not None:
			print(self.gff_dbm._read_data.cache_info(), flush=True)
			self.gff_dbm._read_data.cache_clear()
			bam = BamFile(bamfile_ns)
			self.process_ambiguous_alignments(bam)

		print(self.gff_dbm._read_data.cache_info(), flush=True)
		self.gff_dbm._read_data.cache_clear()
		print(self.gff_dbm._get_tree.cache_info(), flush=True)
		self.gff_dbm._get_tree.cache_clear()

		self.overlap_counter.annotate_counts(bam, self.gff_dbm)
		self.overlap_counter.dump_counts(bam)

		print("Finished.", flush=True)


class AlignmentGroup:
	def __init__(self, aln):
		self.secondaries = list()
		self.primary1, self.primary2 = None, None
		self.qname = aln.qname
		self.add_alignment(aln)
		self.is_ambiguous = aln.mapq == 0

	def add_alignment(self, aln):
		if not aln.is_primary():
			self.secondaries.append(aln)
		elif aln.flag & 0x40:
			self.primary1 = aln
		elif aln.flag & 0x80:
			self.primary2 = aln

	def n_align(self):
		return len(self.secondaries) + int(self.primary1 is not None) + int(self.primary2 is not None)

	def resolve(self, counter, gff_dbm, bam, distmode="all1"):
		hits = dict()
		n_aln = 0
		if self.is_ambiguous:
			unannotated = 0
			alignments = set()
			for aln in [self.primary1, self.primary2] + self.secondaries:
				if aln is not None:
					alignments.add((aln.rid, aln.start, aln.end))
			for rid, start, end in alignments:
				ref = bam.get_reference(rid)[0]
				if gff_dbm.db_index.get(ref) is not None:
					overlaps = gff_dbm.get_overlaps(ref, start, end)
					if overlaps:
						hits.setdefault(rid, set()).update((ovl.begin, ovl.end) for ovl in overlaps)
						n_aln += 1
					else:
						unannotated += 1
				
			counter.update_ambiguous_counts(hits, unannotated, n_aln, gff_dbm, bam, feat_distmode=distmode)
		else:
			# this whole block is for processing unique alignment blocks in a name-sorted bam
			# currently this is not used (-> process_unique_alignments)
			if self.primary1 is not None and self.primary2 is not None and self.primary1.rid == self.primary2.rid:
				start, end = BamFile.calculate_fragment_borders(self.primary1.start, self.primary1.end, self.primary2.start, self.primary2.end)
				alignments = ((start, end),)
			else:
				alignments = (
					(self.primary1.rid, self.primary1.start, self.primary1.end) if self.primary1 else None,
					(self.primary2.rid, self.primary2.start, self.primary2.end) if self.primary2 else None
				)
			for aln in alignments:
				if aln is not None:
					rid, start, end = aln
					ref = bam.get_reference(rid)[0]
					overlaps = gff_dbm.get_overlaps(ref, start, end, cache_data=True)
					counter.update_unique_counts(rid, overlaps)

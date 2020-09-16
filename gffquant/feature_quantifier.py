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

	def annotate_unique_counts(self, bam, gff_dbm):
		'''
		Look up and collate functional annotation for unique hits
		'''
		print("Processing unique counts ...", flush=True)
		t0 = time.time()
		for rid, intervals in self.items():
			ref, _ = bam.get_reference(rid)
			for region, count in intervals.items():
				region_annotation = gff_dbm.get_data(ref, *region)
				for ftype, values in region_annotation: #.items():
					if ftype != "ID":
						for v in values:
							self.featcounts.setdefault(ftype, dict()).setdefault(v, [0, region[1] - region[0] + 1])[0] += count
		self.clear()
		t1 = time.time()
		print("Processed unique counts in {n_seconds}s.".format(n_seconds=t1-t0), flush=True)

	def update_ambiguous_counts(self, hits, gff_dbm):
		for k, v in hits.items():
			self.ambig_counts.setdefault(k, list()).extend(hits)

	@staticmethod
	def calculate_scaling_factor(counts, bam):
		raw_total = sum(counts.values())
		normed_total = sum(count / bam.get_reference(rid)[1] for rid, count in counts.items())
		return raw_total / normed_total
	@staticmethod
	def calculate_feature_scaling_factor(counts):
		raw_total = sum(count for count, flen in counts.values())
		normed_total = sum(count / flen for count, flen in counts.values())
		return raw_total / normed_total

	def dump_counts(self, bam):
		counts_template = "{:d}\t{:.5f}\t{:.5f}"
		seqcount_scaling_factor = OverlapCounter.calculate_scaling_factor(self.seqcounts, bam)

		print("Dumping overlap counters...", flush=True)
		with open("{prefix}.seqname.txt".format(prefix=self.out_prefix), "w") as seq_out:
			for rid, count in self.seqcounts.items():
				seq_id, seq_len = bam.get_reference(rid)
				print(rid, seq_id, seq_len, counts_template.format(*OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)), flush=True, sep="\t", file=seq_out)
		with open("{prefix}.feature_counts.txt".format(prefix=self.out_prefix), "w") as feat_out:
			print("unannotated", self.unannotated_reads, sep="\t", file=feat_out, flush=True)
			for ftype, counts in sorted(self.featcounts.items()):
				print("#{}".format(ftype), file=feat_out, flush=True)
				feature_scaling_factor = OverlapCounter.calculate_feature_scaling_factor(counts)
				for subf, (subf_count, f_len) in sorted(counts.items()):
					print(subf, counts_template.format(*OverlapCounter.normalise_counts(subf_count, f_len, feature_scaling_factor)), flush=True, sep="\t", file=feat_out)

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
		#return IntervalTree.from_tuples(sorted([key[1:] for key in self.loaded_data]))
		return IntervalTree.from_tuples(sorted([key[1:] for key in self._read_data(ref, include_payload=cache_data)]))
	#def _load_data(self, ref, include_payload=False):
	#	if self.loaded_ref != ref:
	#		self.loaded_ref = ref
	#		self.loaded_data = self._read_data(ref, include_payload=include_payload)
	#		self.interval_tree = self._get_tree(ref)  #IntervalTree.from_tuples(sorted([key[1:] for key in self.loaded_data]))

	def get_data(self, ref, start, end):
		#self._load_data(ref, include_payload=True)
		#return self.loaded_data.get((ref, start, end), dict())
		return self._read_data(ref, include_payload=True).get((ref, start, end), dict())

	def get_overlaps(self, ref, start, end, cache_data=False):
		#self._load_data(ref)
		#return self.interval_tree[start:end]
		return self._get_tree(ref, cache_data=cache_data)[start:end]


class FeatureQuantifier:
	def _read_count_config(self, config):
		default = {"multiple": self.default_multiple_alignments, "normalization": self.default_normalization}
		self._count_config = yaml.load(open(config), Loader=yaml.SafeLoader) if config is not None else default

	def __init__(self, gff_db, gff_index, count_config=None, multiple_alignments="unique_only", normalization="scaled"): #, gff_gzipped=True):
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
			qname = aln.qname #get_hash()

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

		self.overlap_counter.annotate_unique_counts(bam, self.gff_dbm)
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

	def resolve(self, counter, gff_dbm, bam):
		hits = Counter()
		if self.is_ambiguous:
			# TODO: scan for duplicates
			for aln in [self.primary1, self.primary2] + self.secondaries:
				if aln is not None:
					ref = bam.get_reference(aln.rid)[0]
					if gff_dbm.db_index.get(ref) is not None:
						overlaps = gff_dbm.get_overlaps(ref, aln.start, aln.end, cache_data=True)
						for ovl in overlaps:
							ann = gff_dbm.get_data(ref, ovl.begin, ovl.end)
							# print(ann, file=sys.stderr, flush=True)
							if ann:
								#try:
								#	del ann["ID"]
								#except:
								#	pass
								#ann_key = tuple(ann.items())
								hits[ann] += 1
								print(ann, hits[ann], file=sys.stderr, flush=True)
								#hits.setdefault(ann_key, set()).add((aln.rid, ovl.begin, ovl.end))
			# counter.update_ambiguous_counts(hits, gff_dbm)
		else:
			if self.primary1 is not None and self.primary2 is not None and self.primary1.rid == self.primary2.rid:
				start, end = BamFile.calculate_fragment_borders(self.primary1.start, self.primary1.end, self.primary2.start, self.primary2.end)
				alignments = ((start, end))
				#overlaps = gff_dbm.get_overlaps(bam.get_reference(self.primary1.rid)[0], start, end, cache_data=True)
				#counter.update_unique_counts(self.primary1.rid, overlaps)
			else:
				alignments = (
					(self.primary1.rid, self.primary1.start, self.primary1.end) if self.primary1 else None,
					(self.primary2.rid, self.primary2.start, self.primary2.end) if self.primary2 else None
				)
			for aln in alignments:
				if aln is not None:
					ref = bam.get_reference(aln.rid)[0]
					overlaps = gff_dbm.get_overlaps(ref, aln.start, aln.end, cache_data=True)
					counter.update_unique_counts(aln.rid, overlaps)

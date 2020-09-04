import os
import time
import sys
import gzip
from collections import Counter
import json
import yaml

from intervaltree import IntervalTree

from bamreader import BamFile

class OverlapCounter(dict):
	@staticmethod
	def normalise_counts(counts, feature_len):
		'''Returns raw, length-normalised, and scaled feature counts.'''
		normalised = counts / feature_len
		return counts, normalised, counts / normalised
	def __init__(self):
		pass
	def update_counts(self, rid, overlaps, n_aln=1):
		hit_list = [(ovl.begin, ovl.end) for ovl in overlaps]
		self.setdefault(rid, list()).append((hit_list, n_aln))
	def dump_counts(self, bam, out_prefix):
		print("Dumping overlap counters...", flush=True)
		with open("{prefix}.seqname.txt".format(prefix=out_prefix), "w") as seq_out:
			for rid, counts in self.items():
				print(rid, *bam.get_reference(rid), *OverlapCounter.normalise_counts(len(counts)), flush=True, sep="\t", file=seq_out)

	"""
	def dump_overlap_counters(self, bam, out_prefix):
		with open(out_prefix + ".feature_counts.txt", "w") as out:
			dump_counter = [dict(), dict()]
			for ref_id in sorted(set(self.unique_counter).union(self.multimap_counter)):
				ref_name = bam.get_reference(ref_id)[0]
				gff_annotation = self._read_gff_data(ref_name, include_payload=True)
				for (start, end) in set(self.unique_counter.get(ref_id, Counter())).union(self.multimap_counter.get(ref_id, Counter())):
					primary_count = self.unique_counter.get(ref_id, Counter())[(start, end)]
					secondary_count = self.multimap_counter.get(ref_id, Counter())[(start, end)]
					for feature_type, values in gff_annotation.get((ref_name, start, end), dict()).items():
						if feature_type != "ID":
							dump_counter[0].setdefault(feature_type, Counter()).update({v: primary_count for v in values})
							dump_counter[1].setdefault(feature_type, Counter()).update({v: secondary_count for v in values})
			json.dump(dump_counter, out, indent="\t")
	"""

class FeatureQuantifier:
	def _read_count_config(self, config):
		default = {"multiple": self.default_multiple_alignments, "normalization": self.default_normalization}
		self._count_config = yaml.load(open(config), Loader=yaml.SafeLoader) if config is not None else default
	def _read_gff_index(self, gff_index):
		self.gff_index = dict()
		for line in open(gff_index, "rt"):
			line = line.strip().split("\t")
			self.gff_index.setdefault(line[0], list()).append(list(map(int, line[1:3])))
	def __init__(self, gff_db, gff_index, count_config=None, multiple_alignments="unique_only", normalization="scaled", gff_gzipped=True):
		self.gff_db = gff_db
		self.default_multiple_alignments = multiple_alignments
		self.default_normalization = normalization
		gff_open = gzip.open if gff_gzipped else open
		self.gff_data = gff_open(gff_db, "rt")

		self._read_gff_index(gff_index)
		self._read_count_config(count_config)
	def _read_gff_data(self, ref_id, include_payload=False):
		gff_annotation = dict()
		for offset, size in self.gff_index.get(ref_id, list()):
			self.gff_data.seek(offset)
			for line in self.gff_data.read(size).strip("\n").split("\n"):
				if not line.startswith("#"):
					line = line.strip().split("\t")
					features = None
					if include_payload:
						features = dict((item.split("=")[0], item.split("=")[1].split(",")) for item in line[8].strip().split(";"))
					key = (line[0], int(line[3]), int(line[4]) + 1)
					gff_annotation[key] = features
		if not gff_annotation:
			print("WARNING: contig {contig} does not have an annotation in the index.".format(contig=ref_id), file=sys.stderr, flush=True)
		return gff_annotation

	def process_unique_cache(self, interval_tree, rid):
		for qname, cached in self.read_cache.items():
			n_aln = len(cached)
			if n_aln > 2:
				print("WARNING: more than two primary alignments for {qname} ({n}). Ignoring.".format(qname=qname, n=n_aln), file=sys.stderr, flush=True)
				continue

			start, end = cached[0][1:] if n_aln == 1 else BamFile.calculate_fragment_borders(*cached[0][1:], *cached[1][1:])
			self.overlap_counter.update_counts(rid, interval_tree[start:end], n_aln=1)

		self.read_cache.clear()

	def update_reference_data(self, ref, rid, current_ref, current_rid, cache, counter, interval_tree):
		current_rid = rid
		if current_ref is not None:
			self.process_unique_cache(interval_tree, rid)
		current_ref = ref
		gff_annotation = self._read_gff_data(current_ref)
		intervals = sorted([key[1:] for key in gff_annotation])
		interval_tree = IntervalTree.from_tuples(intervals)

		return current_rid, current_ref, interval_tree

	def process_bam(self, bamfile, out_prefix):
		bam = BamFile(bamfile)
		multimappers = bam.get_multimappers()
		t0 = time.time()
		self.read_cache = dict()
		self.multimap_cache = dict()
		self.unique_counter, self.multimap_counter = dict(), dict()
		self.overlap_counter = OverlapCounter()

		current_ref, current_rid = None, None
		interval_tree = None

		for aln_count, aln in enumerate(bam.get_alignments(disallowed_flags=0x800)):

			counter = self.unique_counter
			ref = bam.get_reference(aln.rid)[0]
			start, end = aln.start, aln.end

			if ref != current_ref:
				current_rid, current_ref, interval_tree = self.update_reference_data(
					ref,
					aln.rid,
					current_ref,
					current_rid,
					self.read_cache,
					self.unique_counter,
					interval_tree
				)

			mm_count, primaries = multimappers.get(aln.qname, [0, set()])
			if not aln.is_primary() or primaries:
				counter = self.multimap_counter
				if not primaries:
					# ignore secondary alignments without existing primary alignment for now
					print("WARNING: could not find primary alignments for {aln_qname}".format(aln_qname=aln.qname), flush=True, file=sys.stderr)
					continue
				# we need to have the primary alignment flag in the key so we can process the primary alignment properly
				aln_key = (aln.rid, aln.start, aln.end, aln.is_primary())
				# special treatment for the secondary alignments
				# this is stuff i've seen coming out of ngless/bwa
				if not aln_key[-1]:
					mm_count -= 1
					seen_primaries = [k for k in self.multimap_cache.get(aln.qname, set()) if k[-1]]
					if not mm_count and len(seen_primaries) == len(primaries):
						# all primary and secondary alignments of read have been seen
						# TODO process alignments
						self.multimap_cache.pop(aln.qname)
						continue
					multimappers[aln.qname][0] = mm_count

					if aln_key[:-1] in primaries:
						# there are sometimes secondary alignments that seem to be identical to the primary -> ignore
						continue
					if aln_key in self.multimap_cache.get(aln.qname, set()):
						# if r1 and r2 are identical and generate two separate individual alignments -> ignore
						# i don't know why that would happen: very short fragment sequenced twice, i.e. from each direction or r1/r2 mislabeled/duplicated?
						continue
					
				self.multimap_cache.setdefault(aln.qname, set()).add(aln_key)
				continue

			elif aln.is_paired() and aln.rid == aln.rnext:
				# need to keep track of read pairs to avoid dual counts
				# check if the mate has already been seen
				mates = self.read_cache.setdefault(aln.qname, list())
				if mates:
					# if it has, calculate the total fragment size and remove the pair from the cache
					if aln.rnext != mates[0][0]:
						print("WARNING: alignment {qname} seems to be corrupted: {aln1} {aln2}".format(qname=aln.qname, aln1=str(aln), aln2=str(mates[0])), flush=True, file=sys.stderr)
						continue # i don't think this ever happens
					else:
						start, end = BamFile.calculate_fragment_borders(aln.start, aln.end, *mates[0][1:])
					self.read_cache.pop(aln.qname)
				else:
					# otherwise cache the first encountered mate and advance to the next read
					mates.append((aln.rid, aln.start, aln.end))
					continue

			# at this point only single-end reads and merged pairs should be processed here
			self.overlap_counter.update_counts(aln.rid, interval_tree[start:end], n_aln=1)

		# clear cache
		self.process_unique_cache(interval_tree, current_rid)
		t1 = time.time()

		print("Processed {n_align} primary alignments in {n_seconds:.3f}s.".format(
			n_align=aln_count, n_seconds=t1-t0), flush=True
		)

		self.overlap_counter.dump_counts(bam, out_prefix)

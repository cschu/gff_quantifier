import os
import time
import sys
import gzip
from collections import Counter
import json
import yaml
from datetime import datetime

from intervaltree import IntervalTree

from bamreader import BamFile

class OverlapCounter(dict):
	@staticmethod
	def normalise_counts(counts, feature_len):
		'''Returns raw, length-normalised, and scaled feature counts.'''
		normalised = counts / feature_len
		scaled = normalised * counts / normalised #? -> that's what it seems to do in ngless ???
		fpkm = ((normalised * 10e9) / counts) * normalised  # ?
		return counts, normalised, scaled, fpkm
	def __init__(self, out_prefix):
		self.out_prefix = out_prefix
		self.seqcounts = Counter()
		self.featcounts = dict()
		pass
	def update_counts(self, rid, overlaps, n_aln=1):
		'''
		Generates a hit list from the overlaps resulting from an intervaltree query,
		adds the number of alternative alignments and stores the results for each reference sequence.
		'''
		hit_list = [(ovl.begin, ovl.end) for ovl in overlaps]
		self.setdefault(rid, list()).append((hit_list, n_aln))

	def process_counts(self, bam, gff_source):
		print("Processing counts ...", flush=True)
		t0 = time.time()
		# first pass: process unique mappers
		for rid, counts in self.items():
			ref, reflen = bam.get_reference(rid)
			gff_annotation = gff_source.read_gff_data(ref, include_payload=True)
			mcounts = list()
			for hits, n_aln in counts:
				if n_aln == 1:
					# we don't care about hits to overlapping features when counting ref-hits
					self.seqcounts[rid] += 1
				if n_aln == 1 and len(hits) == 1:
					# add single count to each subfeature
					for ftype, values in gff_annotation.get((ref, *hits[0]), dict()).items():
						if ftype != "ID":
							for v in values:
								self.featcounts.setdefault(ftype, dict()).setdefault(v, [0, hits[0][1]-hits[0][0]+1])[0] += 1
				else:
					mcounts.append((hits, n_aln))
			self[rid] = mcounts
		t1 = time.time()
		print("Processed unique hits in {n_seconds}s.".format(n_seconds=t1-t0), flush=True)
	def dump_counts(self, bam):
		print("Dumping overlap counters...", flush=True)
		with open("{prefix}.seqname.txt".format(prefix=self.out_prefix), "w") as seq_out:
			for rid, count in self.seqcounts.items():
				seq_id, seq_len = bam.get_reference(rid)
				print(rid, seq_id, seq_len, "{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}".format(*OverlapCounter.normalise_counts(count, seq_len)), flush=True, sep="\t", file=seq_out)
		with open("{prefix}.feature_counts.txt".format(prefix=self.out_prefix), "w") as feat_out:
			for ftype, counts in sorted(self.featcounts.items()):
				print("#{}".format(ftype), file=feat_out, flush=True)
				for subf, (subf_count, f_len) in sorted(counts.items()):
					print(subf, "{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}".format(*OverlapCounter.normalise_counts(subf_count, f_len)), flush=True, sep="\t", file=feat_out)


class FeatureQuantifier:
	def _read_count_config(self, config):
		default = {"multiple": self.default_multiple_alignments, "normalization": self.default_normalization}
		self._count_config = yaml.load(open(config), Loader=yaml.SafeLoader) if config is not None else default
	def read_gff_index(self, gff_index):
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

		self.read_gff_index(gff_index)
		self._read_count_config(count_config)
	def read_gff_data(self, ref_id, include_payload=False):
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
		if not gff_annotation and not include_payload:
			print("WARNING: contig {contig} does not have an annotation in the index.".format(contig=ref_id), file=sys.stderr, flush=True)
		return gff_annotation

	def process_unique_cache(self, interval_tree, rid):
		for qname, cached in self.umap_cache.items():
			n_aln = len(cached)
			if n_aln > 2:
				print("WARNING: more than two primary alignments for {qname} ({n}). Ignoring.".format(qname=qname, n=n_aln), file=sys.stderr, flush=True)
				continue

			start, end = cached[0][1:] if n_aln == 1 else BamFile.calculate_fragment_borders(*cached[0][1:], *cached[1][1:])
			self.overlap_counter.update_counts(rid, interval_tree[start:end], n_aln=1)

		self.umap_cache.clear()

	def update_reference_data(self, ref, rid, current_ref, interval_tree):


		if current_ref is not None:
			self.process_unique_cache(interval_tree, rid)
		current_ref = ref
		gff_annotation = self.read_gff_data(current_ref)
		intervals = sorted([key[1:] for key in gff_annotation])
		interval_tree = IntervalTree.from_tuples(intervals)

		return rid, current_ref, interval_tree

	def process_bam(self, bamfile, out_prefix):
		bam = BamFile(bamfile)
		# TODO multimappers
		# multimappers = bam.get_multimappers()
		multimappers = dict()
		t0 = time.time()
		self.umap_cache, self.mmap_cache = dict(), dict()
		self.overlap_counter = OverlapCounter(out_prefix)

		current_ref, current_rid = None, None
		interval_tree = None

		for aln_count, aln in bam.get_alignments(disallowed_flags=0x800):

			ref = bam.get_reference(aln.rid)[0]
			start, end = aln.start, aln.end
			qname = aln.get_hash()

			if ref != current_ref:
				current_rid, current_ref, interval_tree = self.update_reference_data(
					ref, aln.rid, current_ref, interval_tree)
				print("{time}\tNew reference: {ref} ({rid}/{n_ref}). {n_aln} alignments processed.".format(
					time=datetime.now().strftime("%m/%d/%Y,%H:%M:%S"),
					ref=current_ref,
					rid=current_rid,
					n_ref=bam.n_references(),
					n_aln=aln_count),
					file=sys.stderr, flush=True)


			# if not aln.is_primary() or aln.mapq == 0:
			if aln.mapq == 0:
				continue # TODO: re-enable multimappers
				mm_count, primaries = multimappers.get(qname, [0, set()])
				if not primaries:
					# ignore secondary alignments without existing primary alignment for now
					#print("WARNING: could not find primary alignments for {aln_qname}".format(aln_qname=qname), flush=True, file=sys.stderr)
					continue
				# we need to have the primary alignment flag in the key so we can process the primary alignment properly
				aln_key = (aln.rid, aln.start, aln.end, aln.is_primary())
				# special treatment for the secondary alignments
				# this is stuff i've seen coming out of ngless/bwa
				if not aln_key[-1]:
					mm_count -= 1
					seen_primaries = [k for k in self.mmap_cache.get(qname, set()) if k[-1]]
					if not mm_count and len(seen_primaries) == len(primaries):
						# all primary and secondary alignments of read have been seen
						# TODO process alignments
						self.mmap_cache.pop(qname)
						continue
					multimappers[qname][0] = mm_count

					if aln_key[:-1] in primaries:
						# there are sometimes secondary alignments that seem to be identical to the primary -> ignore
						continue
					if aln_key in self.mmap_cache.get(qname, set()):
						# if r1 and r2 are identical and generate two separate individual alignments -> ignore
						# i don't know why that would happen: very short fragment sequenced twice, i.e. from each direction or r1/r2 mislabeled/duplicated?
						continue
					
				self.mmap_cache.setdefault(qname, set()).add(aln_key)
				continue

			elif aln.is_paired() and aln.rid == aln.rnext:
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
					# self.umap_cache.pop(qname)
					del self.umap_cache[qname]
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

		self.overlap_counter.process_counts(bam, self)
		self.overlap_counter.dump_counts(bam)
		print("Finished.", flush=True)

import os
import time
import sys
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
	def update_unique_counts(self, rid, overlaps):
		'''
		Generates a hit list from the overlaps resulting from an intervaltree query,
		adds the number of alternative alignments and stores the results for each reference sequence.
		'''
		self.setdefault(rid, Counter()).update((ovl.begin, ovl.end) for ovl in overlaps)
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
				for ftype, values in region_annotation.items():
					if ftype != "ID":
						for v in values:
							self.featcounts.setdefault(ftype, dict()).setdefault(v, [0, region[1] - region[0] + 1])[0] += count
			del self[rid]
		self.clear()
		t1 = time.time()
		print("Processed unique counts in {n_seconds}s.".format(n_seconds=t1-t0), flush=True)
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
	def _read_data(self, ref_id, include_payload=False):
		gff_annotation = dict()
		for offset, size in self.db_index.get(ref_id, list()):
			self.db.seek(offset)
			for line in self.db.read(size).strip("\n").split("\n"):
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
	def _load_data(self, ref, include_payload=False):
		if self.loaded_ref != ref:
			self.loaded_ref = ref
			self.loaded_data = self._read_data(ref, include_payload=include_payload)
			self.interval_tree = IntervalTree.from_tuples(sorted([key[1:] for key in self.loaded_data]))

	def get_data(self, ref, start, end):
		self._load_data(ref, include_payload=True)
		return self.loaded_data.get((ref, start, end), dict())

	def get_overlaps(self, ref, start, end):
		self._load_data(ref)
		return self.interval_tree[start:end]


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

	def process_unique_alignments(self, bam):
		t0 = time.time()
		current_ref, current_rid = None, None

		for aln_count, aln in bam.get_alignments(allow_multiple=False, disallowed_flags=0x800):
			ref = bam.get_reference(aln.rid)[0]
			start, end = aln.start, aln.end
			qname = aln.get_hash()

			if ref != current_ref:
				print("{time}\tNew reference: {ref} ({rid}/{n_ref}). {n_aln} alignments processed.".format(
					time=datetime.now().strftime("%m/%d/%Y,%H:%M:%S"),
					ref=ref, rid=aln.rid, n_ref=bam.n_references, n_aln=aln_count),
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

	def process_bam(self, bamfile, out_prefix):
		self.overlap_counter = OverlapCounter(out_prefix)
		bam = BamFile(bamfile)

		self.process_unique_alignments(bam)
		#bam.rewind()
		#self.process_multiple_alignments(bam)

		self.overlap_counter.annotate_unique_counts(bam, self.gff_dbm)
		self.overlap_counter.dump_counts(bam)
		print("Finished.", flush=True)

	def process_multiple_alignments(self, bam):
		multiples = list() #dict()
		size, n_aln = 0, 0
		t0 = time.time()
		read_lut = dict()
		for aln_count, aln in bam.get_alignments(allow_unique=False, disallowed_flags=0x800):
			qname, v = aln.get_hash(), (aln.rid, aln.start, aln.end, aln.is_primary())
			size += sys.getsizeof(qname) + sum(map(sys.getsizeof, v))
			n_aln += 1
			read_id = read_lut.setdefault(qname, len(read_lut))
			#multiples.setdefault(aln.get_hash(), list()).append((aln.rid, aln.start, aln.end, aln.is_primary()))
			multiples.append((read_id, aln.rid, aln.start, aln.end, aln.is_primary()))
		size += sys.getsizeof(multiples)
		t1 = time.time()
		print("Processed {n_align} secondary alignments (size={size}b) in {n_seconds:.3f}s.".format(
			n_align=n_aln, n_seconds=t1-t0, size=size), flush=True)

		n_no_primary = 0
		#for alignments in multiples: #.values():
		#	if all(not is_primary for _, _, _, is_primary in alignments):
		#		n_no_primary += 1
		#print("{n_no_primary} secondary alignments don't have a primary alignment.".format(n_no_primary=n_no_primary), flush=True)

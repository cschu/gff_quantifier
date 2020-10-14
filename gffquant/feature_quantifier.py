import csv
import os
import time
import sys
from collections import Counter
from datetime import datetime

from gffquant.bamreader import BamFile
from gffquant.gff_dbm import GffDatabaseManager
from gffquant.overlap_counter import OverlapCounter

SUPPL_ALN_FLAG = 0x800
WEIRD_MASK_FLAG = 0x1c0 # this masks everything except first/second in pair and secondary aln
MATE_UNMAPPED_FLAG = 0x8
FIRST_IN_PAIR_FLAG = 0x40
SECOND_IN_PAIR_FLAG = 0x80

class FeatureQuantifier:

	def __init__(self, gff_db, gff_index, out_prefix="gffquant", ambig_mode="unique_only"):
		self.gff_dbm = GffDatabaseManager(gff_db, gff_index)
		self.umap_cache = dict()
		self.overlap_counter = OverlapCounter(out_prefix)
		self.out_prefix = out_prefix
		self.ambig_mode = ambig_mode

	def process_unique_cache(self, rid, ref):
		for qname, uniq_aln in self.umap_cache.items():
			n_aln = len(uniq_aln)
			if n_aln > 2:
				print("WARNING: more than two primary alignments for {qname} ({n}). Ignoring.".format(qname=qname, n=n_aln), file=sys.stderr, flush=True)
				continue

			start, end = uniq_aln[0][1:] if n_aln == 1 else BamFile.calculate_fragment_borders(*uniq_aln[0][1:], *uniq_aln[1][1:])
			self.overlap_counter.update_unique_counts(rid, self.gff_dbm.get_overlaps(ref, start, end))

		self.umap_cache.clear()

	def process_alignments(self, bam):
		t0 = time.time()
		current_ref, current_rid = None, None

		with open(self.out_prefix + ".ambig_tmp.txt", "wt") as ambig_out_tmp:
			for aln_count, aln in bam.get_alignments(allow_multiple=self.ambig_mode != "unique_only", allow_unique=True, disallowed_flags=SUPPL_ALN_FLAG):
				start, end = aln.start, aln.end
				qname = aln.qname

				if aln.rid != current_rid:
					if current_rid is not None:
						self.process_unique_cache(current_rid, current_ref)
					current_ref, current_rid = bam.get_reference(aln.rid)[0], aln.rid
					
					print("{time}\tNew reference: {ref} ({rid}/{n_ref}). {n_aln} alignments processed.".format(
						time=datetime.now().strftime("%m/%d/%Y,%H:%M:%S"),
						ref=current_ref, rid=aln.rid, n_ref=bam.n_references(), n_aln=aln_count),
						file=sys.stderr, flush=True)

				if aln.is_ambiguous() and self.ambig_mode == "1overN":
					overlaps = self.gff_dbm.get_overlaps(current_ref, start, end)
					if not overlaps:
						print(aln.qname, aln_count, -1, -1, -1, aln.flag, file=ambig_out_tmp, sep="\t")
					for ovl in overlaps:
						print(aln.qname, aln_count, aln.rid, ovl.begin, ovl.end, aln.flag, file=ambig_out_tmp, sep="\t")
				elif aln.is_unique() and aln.is_paired() and aln.rid == aln.rnext and not aln.flag & MATE_UNMAPPED_FLAG:
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

				# at this point only single-end reads and merged pairs should be processed here ( + ambiguous reads in "all1" mode )
				self.overlap_counter.update_unique_counts(aln.rid, self.gff_dbm.get_overlaps(current_ref, start, end))

		self.process_unique_cache(current_rid, current_ref)
		t1 = time.time()

		print("Processed {n_align} alignments in {n_seconds:.3f}s.".format(
		    n_align=aln_count, n_seconds=t1-t0), flush=True
		)

	def _read_ambiguous_alignments(self, ambig_in):
		read_ids, ambig_aln = dict(), list()
		for qname, *read_data in csv.reader(ambig_in, delimiter="\t"):
			qname = read_ids.setdefault(qname, len(read_ids))
			ambig_aln.append((qname, ) + tuple(map(int, read_data)))

		return ambig_aln

	def _process_ambiguous_aln_groups(self, ambig_aln, bam):
		n_align, current_group = 0, None
		for aln in sorted(ambig_aln):
			if current_group is None or current_group.qname != aln[0]:
				if current_group:
					n_align += current_group.n_align()
					current_group.resolve(self.overlap_counter, self.gff_dbm, bam)
				current_group = AmbiguousAlignmentGroup(aln)
			else:
				current_group.add_alignment(aln)

		if current_group is not None:
			n_align += current_group.n_align()
			current_group.resolve(self.overlap_counter, self.gff_dbm, bam)

		return n_align
		

	def process_data(self, bamfile):
		bam = BamFile(bamfile)
		self.process_alignments(bam)

		print(self.gff_dbm._read_data.cache_info(), flush=True)
		self.gff_dbm._read_data.cache_clear()
		print(self.gff_dbm._get_tree.cache_info(), flush=True)
		self.gff_dbm._get_tree.cache_clear()

		if self.ambig_mode != "unique_only":
			# b'ST-K00119:83:HH22HBBXX:1:2203:9912:18107\x00' 4       126964  129574  2       0       256
			t0 = time.time()
			with open(self.out_prefix + ".ambig_tmp.txt") as ambig_in:
				ambig_aln = self._read_ambiguous_alignments(ambig_in)

			n_align = self._process_ambiguous_aln_groups(ambig_aln, bam)
			os.remove(self.out_prefix + ".ambig_tmp.txt")

			t1 = time.time()
			print("Processed {n_align} secondary alignments in {n_seconds:.3f}s.".format(
				n_align=n_align, n_seconds=t1-t0), flush=True)

		self.overlap_counter.annotate_counts(bam, self.gff_dbm)
		self.overlap_counter.dump_counts(bam)

		print("Finished.", flush=True)


class AmbiguousAlignmentGroup:
	def __init__(self, aln):
		self.secondaries = list()
		self.primary1, self.primary2 = None, None
		self.qname = aln[0]
		self.uniq_read_ids = set()
		self.unannotated = 0
		self.add_alignment(aln)

	# (read, int(aln), int(rid), int(start), int(end), int(flag))
	# int_qname        8       -1      -1      -1      256
	def add_alignment(self, aln):
		short_aln = aln[1:] #(aln.rid, aln.start, aln.end)
		if aln[0] == -1:
			self.unannotated += 1
		elif aln[-1] & FIRST_IN_PAIR_FLAG and self.primary1 is None:
			self.primary1 = short_aln
		elif aln[-1] & SECOND_IN_PAIR_FLAG and self.primary2 is None:
			self.primary2 = short_aln
		else:
			self.secondaries.append(short_aln)
		self.uniq_read_ids.add(aln[0])

	def n_align(self):
		return len(self.uniq_read_ids.difference({-1}))

	def resolve(self, counter, gff_dbm, bam, distmode="all1"):
		hits = dict()

		alignments = {aln[1:] for aln in [self.primary1, self.primary2] + self.secondaries if aln is not None and aln[1] != -1}
		for rid, start, end, flag in alignments:
			hits.setdefault(rid, set()).add((start, end))

		counter.update_ambiguous_counts(hits, self.unannotated, self.n_align(), gff_dbm, bam, feat_distmode=distmode)

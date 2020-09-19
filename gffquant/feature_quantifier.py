import os
import time
import sys
from collections import Counter
from datetime import datetime

from gffquant.bamreader import BamFile
from gffquant.gff_dbm import GffDatabaseManager
from gffquant.overlap_counter import OverlapCounter


class FeatureQuantifier:

	def __init__(self, gff_db, gff_index, out_prefix="gffquant"):
		self.gff_dbm = GffDatabaseManager(gff_db, gff_index)
		self.umap_cache = dict()
		self.overlap_counter = OverlapCounter(out_prefix)

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

	def process_data(self, bamfile, bamfile_ns=None):
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

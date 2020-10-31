import csv
import os
import time
import sys
from collections import Counter
from datetime import datetime
import contextlib

import numpy
import pandas

from gffquant.bamreader import BamFile, SamFlags
from gffquant.gff_dbm import GffDatabaseManager
from gffquant.overlap_counter import OverlapCounter

SUPPL_ALN_FLAG = 0x800 # we don't allow supplementary alignments
MATE_UNMAPPED_FLAG = 0x8
FIRST_IN_PAIR_FLAG = 0x40
SECOND_IN_PAIR_FLAG = 0x80
REVCOMP_ALIGNMENT_FLAG = 0x10

DEBUG = True

class AmbiguousAlignmentRecordKeeper:
	"""
	This class takes care of the specific record keeping for ambiguous alignments.
	This includes:
		- counting of annotated / unannotated reads
		- assignment of integer read ids to save space
		- writing of annotated read information to a temporary file
		(this last one is a painful workaround, the alternative would be to
		 split the bam file into unique and ambiguous alignments, with the latter sorted by name)
		 The task here is to save memory, which now comes at the cost of a bit of (temporary) disk space.

	The class is used as a contextmanager, it is only instanced when ambiguous alignments require special treatment.
	"""

	def __init__(self, prefix):
		self.annotated = set()
		self.unannotated = set()
		self.readids = dict()
		self.dumpfile = prefix + ".ambig_tmp.txt"
		self.ambig_dump = open(self.dumpfile, "at")
	def __enter__(self):
		return self
	def __exit__(self, exc_type, exc_val, exc_tb):
		self.ambig_dump.close()
	def process_alignment(self, aln, aln_count, overlaps=None):
		"""
		- obtains/assigns a unique integer (correlating to the read id) to an alignment for alignment group identification
		- updates the unannotated/annotated records, depending on whether the alignment could be annotated
		- writes the relevant alignment information to disk if alignment could be annotated
		"""
		qname_id = self.readids.setdefault(aln.qname, len(self.readids))
		if not overlaps:
			self.unannotated.add(qname_id)
		else:
			self.annotated.add(qname_id)
			for ovl in overlaps:
				print(qname_id, aln_count, aln.rid, ovl.begin, ovl.end, aln.flag, file=self.ambig_dump, sep="\t")
	def n_unannotated(self):
		""" returns the number of unannotated reads (all reads that didn't align to any annotated region) """
		return len(self.unannotated.difference(self.annotated))



class FeatureQuantifier:
	TRUE_AMBIG_MODES = ("dist1", "1overN")

	def __init__(self, gff_db, gff_index, out_prefix="gffquant", ambig_mode="unique_only"):
		self.gff_dbm = GffDatabaseManager(gff_db, gff_index)
		self.umap_cache = dict()
		self.overlap_counter = OverlapCounter(out_prefix)
		self.out_prefix = out_prefix
		self.ambig_mode = ambig_mode
		print("Ambig mode:", self.ambig_mode)

	def allow_ambiguous_alignments(self):
		return self.ambig_mode != "unique_only"
	def treat_ambiguous_as_unique(self):
		return self.ambig_mode not in FeatureQuantifier.TRUE_AMBIG_MODES
	def require_ambig_bookkeeping(self):
		return self.allow_ambiguous_alignments() and not self.treat_ambiguous_as_unique()

	def process_unique_cache(self, rid, ref):
		for qname, uniq_aln in self.umap_cache.items():
			n_aln = len(uniq_aln)
			if n_aln == 1:
				start, end, flag = uniq_aln[0][1:]
				rev_strand = SamFlags.is_reverse_strand(flag) # & REVCOMP_ALIGNMENT_FLAG
			elif n_aln == 2:
				start, end = BamFile.calculate_fragment_borders(*uniq_aln[0][1:-1], *uniq_aln[1][1:-1])
				rev_strand = None #TODO: add strand-specific handling by RNAseq protocol
			else:
				print("WARNING: more than two primary alignments for {qname} ({n}). Ignoring.".format(qname=qname, n=n_aln), file=sys.stderr, flush=True)
				continue

			self.overlap_counter.update_unique_counts(rid, self.gff_dbm.get_overlaps(ref, start, end), rev_strand=rev_strand)

		self.umap_cache.clear()

	def process_alignments(self, bam):
		"""
		Reads from a position-sorted bam file are processed in batches according to the reference sequence they were aligned to.
		This allows a partitioning of the reference annotation (which saves memory and time in case of large reference data sets).
		Ambiguous reads are dumped to disk for dist1 and 1overN distribution modes and processed in a follow-up step.
		"""
		t0 = time.time()
		current_ref, current_rid = None, None

		bam_stream = bam.get_alignments(
			allow_multiple=self.allow_ambiguous_alignments(),
			allow_unique=True, disallowed_flags=SUPPL_ALN_FLAG)

		ambig_bookkeeper = AmbiguousAlignmentRecordKeeper(self.out_prefix) if self.require_ambig_bookkeeping() else contextlib.nullcontext()
		with ambig_bookkeeper:
			for aln_count, aln in bam_stream:
				start, end, rev_strand = aln.start, aln.end, aln.is_reverse() #flag & REVCOMP_ALIGNMENT_FLAG

				if aln.rid != current_rid:
					# if a new reference sequence is encountered, empty the alignment cache (if needed)
					if current_rid is not None:
						self.process_unique_cache(current_rid, current_ref)
					current_ref, current_rid = bam.get_reference(aln.rid)[0], aln.rid

					print("{time}\tNew reference: {ref} ({rid}/{n_ref}). {n_aln} alignments processed.".format(
						time=datetime.now().strftime("%m/%d/%Y,%H:%M:%S"),
						ref=current_ref, rid=aln.rid, n_ref=bam.n_references(), n_aln=aln_count),
						file=sys.stderr, flush=True)

				if aln.is_ambiguous() and self.require_ambig_bookkeeping():
					# if ambiguous alignments are not treated as individual alignments (dist1 or 1overN mode)
					# then the alignment processing is deferred to the bookkeeper
					overlaps = self.gff_dbm.get_overlaps(current_ref, start, end)
					ambig_bookkeeper.process_alignment(aln, aln_count, overlaps=overlaps)
					start, end, rev_strand = None, None, None
				elif aln.is_unique() and aln.is_paired() and aln.rid == aln.rnext and not aln.flag & MATE_UNMAPPED_FLAG:
					# if a read belongs to a properly-paired (both mates align to the same reference sequence), unique alignment
					# it can only be processed if the mate has already been encountered, otherwise it is stored
					start, end, rev_strand = self._process_unique_properly_paired_alignment(aln)

				# at this point only single-end reads, 'improper' and merged pairs should be processed here ( + ambiguous reads in "all1" mode )
				if start is not None:
					self.overlap_counter.update_unique_counts(aln.rid, self.gff_dbm.get_overlaps(current_ref, start, end), rev_strand=rev_strand)

			self.process_unique_cache(current_rid, current_ref)
			t1 = time.time()

			print("Processed {n_align} alignments in {n_seconds:.3f}s.".format(
				n_align=aln_count, n_seconds=t1-t0), flush=True
			)

		if self.require_ambig_bookkeeping():
			return ambig_bookkeeper.n_unannotated(), ambig_bookkeeper.dumpfile 

		return 0, None

	def _process_unique_properly_paired_alignment(self, aln):
		# need to keep track of read pairs to avoid dual counts
		# check if the mate has already been seen
		mates = self.umap_cache.setdefault(aln.qname, list())
		start, end, rev_strand = None, None, None #TODO this requires some extra magic for strand-specific paired-end RNAseq
		if mates:
			# if it has, calculate the total fragment size and remove the pair from the cache
			if aln.rnext != mates[0][0]:
				print("WARNING: alignment {qname} seems to be corrupted: {aln1} {aln2}".format(qname=aln.qname, aln1=str(aln), aln2=str(mates[0])), flush=True, file=sys.stderr)
				# continue # i don't think this ever happens
			else:
				# for a proper pair, calculate the coordinates of the whole fragment
				#TODO could use a length cutoff here (treat mates that align too far apart as individual alignments)
				start, end = BamFile.calculate_fragment_borders(aln.start, aln.end, *mates[0][1:-1])
				rev_strand = None
				del self.umap_cache[aln.qname]
		else:
			# otherwise cache the first encountered mate and advance to the next read
			mates.append((aln.rid, aln.start, aln.end, aln.flag))

		return start, end, rev_strand

	def _read_ambiguous_alignments(self, ambig_in):
		""" reads the dumped ambig alignments and returns them sorted by group """
		# using a pandas dataframe/numpy array saves us a lot of memory
		ambig_aln = pandas.read_csv(ambig_in, sep="\t", header=None)
		return ambig_aln.sort_values(axis=0, by=0)

	def _process_ambiguous_aln_groups(self, ambig_aln, bam):
		""" iterates over name-sorted (grouped) alignments and processes groups individually """
		n_align, current_group = 0, None
		# losing some speed due to the row-iterations here
		# but ok for now
		for aln in ambig_aln.itertuples(index=False, name=None):
			if current_group is None or current_group.qname != aln[0]:
				if current_group:
					n_align += current_group.n_align()
					current_group.resolve(self.overlap_counter, self.gff_dbm, bam, distmode=self.ambig_mode)
				current_group = AmbiguousAlignmentGroup(aln)
			else:
				current_group.add_alignment(aln)

		if current_group is not None:
			n_align += current_group.n_align()
			current_group.resolve(self.overlap_counter, self.gff_dbm, bam, distmode=self.ambig_mode)

		return n_align

	def process_data(self, bamfile, strand_specific=False):
		""" processes one position-sorted bamfile """
		bam = BamFile(bamfile)
		# first pass: process uniqs and dump ambigs (if required)
		unannotated_ambig, ambig_dumpfile = self.process_alignments(bam)
		self.gff_dbm.clear_caches()

		# second pass: process ambiguous alignment groups
		if self.ambig_mode in FeatureQuantifier.TRUE_AMBIG_MODES:
			t0 = time.time()

			ambig_aln = self._read_ambiguous_alignments(ambig_dumpfile)
			n_align = self._process_ambiguous_aln_groups(ambig_aln, bam)

			if not DEBUG:
				os.remove(ambig_dumpfile)

			t1 = time.time()
			print("Processed {n_align} secondary alignments in {n_seconds:.3f}s.".format(
				n_align=n_align, n_seconds=t1-t0), flush=True)

		self.overlap_counter.annotate_counts(bam, self.gff_dbm, strand_specific=strand_specific)
		self.overlap_counter.unannotated_reads += unannotated_ambig
		self.overlap_counter.dump_counts(bam, strand_specific=strand_specific)

		print("Finished.", flush=True)


class AmbiguousAlignmentGroup:

	"""
	Represents a group of ambiguous alignments of a single read/read pair.
	bwa -A assigns a primary read (pair) and flags all others as secondary alignments (0x100)
	due to the ngless filtering, we also have the case that the primary was filtered out

	Caveat: paired-end information is currently not considered for ambiguous alignment groups.
	This means that overlapping mates will result in duplicate counts.
	"""

	def __init__(self, aln):
		self.secondaries = list() # these are the secondary alignments
		self.primary1, self.primary2 = None, None
		self.qname = aln[0]
		self.uniq_alignments = set() # this is the set of alignments that can be annotated
		self.unannotated = 0
		self.add_alignment(aln)

	def add_alignment(self, aln):
		qname, flag = aln[0], aln[-1]
		short_aln = tuple(aln[2:])
		if short_aln[0] == -1:
			self.unannotated += 1
		elif flag & FIRST_IN_PAIR_FLAG and self.primary1 is None:
			self.primary1 = short_aln
		elif flag & SECOND_IN_PAIR_FLAG and self.primary2 is None:
			self.primary2 = short_aln
		else:
			self.secondaries.append(short_aln)
		if short_aln[0] != -1:
			self.uniq_alignments.add(short_aln[:-1])

	def n_align(self):
		return len(self.uniq_alignments.union((self.primary1, self.primary2)).difference({None}))

	def resolve(self, counter, gff_dbm, bam, distmode="all1"):

		#print("RESOLVING", self.qname, self.primary1, self.primary2, self.unannotated, self.n_align(), list(self.uniq_alignments)[:10])

		hits = dict()

		if self.primary1 is not None and self.primary2 is not None and self.primary1[:-1] == self.primary2[:-1]:
			self.primary2 = None

		alignments = set([self.primary1, self.primary2]).union(self.secondaries).difference({None})
		for rid, start, end, flag in alignments:
			# print(rid, bam.get_reference(rid), start, end, flag)
			hits.setdefault(rid, set()).add((start, end, SamFlags.is_reverse_strand(flag))) # & REVCOMP_ALIGNMENT_FLAG))

		counter.update_ambiguous_counts(hits, self.n_align(), self.unannotated, gff_dbm, bam, feat_distmode=distmode)

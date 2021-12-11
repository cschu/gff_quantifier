import csv
import os
import time
import sys
from datetime import datetime
import contextlib

import pandas

from gffquant.bamreader import BamFile, SamFlags
from gffquant.gff_dbm import GffDatabaseManager
# from gffquant.overlap_counter import OverlapCounter
from gffquant.alignment import AmbiguousAlignmentRecordKeeper, AmbiguousAlignmentGroup, PairedEndAlignmentCache

from gffquant.counters import CountManager
 

class FeatureQuantifier:
	TRUE_AMBIG_MODES = ("dist1", "1overN")

	@staticmethod
	def read_ambiguous_alignments(ambig_in):
		""" reads the dumped ambig alignments and returns them sorted by group """
		# using a pandas dataframe/numpy array saves us a lot of memory
		ambig_aln = pandas.read_csv(ambig_in, sep="\t", header=None)
		return ambig_aln.sort_values(axis=0, by=0)

	def __init__(
		self, db=None, db_index=None, out_prefix="gffquant", ambig_mode="unique_only",
		reference_type="genome", strand_specific=False, emapper_version=None, debugmode=False
	):
		self.gff_dbm = GffDatabaseManager(db, reference_type, db_index=db_index, emapper_version=emapper_version)
		self.umap_cache = PairedEndAlignmentCache()
		self.ambig_cache = PairedEndAlignmentCache(ambig_alignments=True)
		self.do_overlap_detection = reference_type in ("genome", "domain")
		#self.overlap_counter = OverlapCounter(
		#	out_prefix, self.gff_dbm,
		#	do_overlap_detection=self.do_overlap_detection, strand_specific=strand_specific,
		#	feature_distribution=ambig_mode
		#)
		self.count_manager = CountManager(
			distribution_mode=ambig_mode,
			region_counts = reference_type in ("genome", "domain"),
			strandedness_required=strand_specific and reference_type not in ("genome", "domain")
		)
		self.out_prefix = out_prefix
		self.ambig_mode = ambig_mode
		self.bamfile = None
		self.debugmode = debugmode
		print("Ambig mode:", self.ambig_mode, flush=True)

	def allow_ambiguous_alignments(self):
		return self.ambig_mode != "unique_only"

	def treat_ambiguous_as_unique(self):
		return self.ambig_mode not in FeatureQuantifier.TRUE_AMBIG_MODES

	def require_ambig_bookkeeping(self):
		return self.allow_ambiguous_alignments() and not self.treat_ambiguous_as_unique()

	def process_alignments_sameref(self, ref, alignments):
		for rid, start, end, rev_strand in alignments:
			if self.do_overlap_detection:
				overlaps, coverage = self.gff_dbm.get_overlaps(ref, start, end)
				hits = {
					(ovl.begin, ovl.end, rev_strand, cstart, cend)
					for ovl, (cstart, cend) in zip(overlaps, coverage)
				}
				aln_count = int(bool(overlaps))
			else:
				# hits = {(-1, -1, rev_strand, None, None)}
				hits = {(None, None, rev_strand, None, None)}
				aln_count = 1

			yield ({rid: hits}, aln_count, 1 - aln_count)

	def process_caches(self, ref):
		self.count_manager.update_unique_counts(
			self.process_alignments_sameref(ref, self.umap_cache.empty_cache())
		)
		self.count_manager.update_ambiguous_counts(
			self.process_alignments_sameref(ref, self.ambig_cache.empty_cache())
		)

	def process_alignments(self, min_identity=None, min_seqlen=None):
		"""
		Reads from a position-sorted bam file are processed in batches according to the reference
		sequence they were aligned to. This allows a partitioning of the reference annotation
		(which saves memory and time in case of large reference data sets).
		Ambiguous reads are dumped to disk for dist1 and 1overN distribution modes and processed
		in a follow-up step.
		"""
		t0 = time.time()
		current_ref, current_rid = None, None

		bam_stream = self.bamfile.get_alignments(
			allow_multiple=self.allow_ambiguous_alignments(),
			allow_unique=True,
			disallowed_flags=SamFlags.SUPPLEMENTARY_ALIGNMENT,
			min_identity=min_identity,
			min_seqlen=min_seqlen
		)

		if self.require_ambig_bookkeeping():
			ambig_bookkeeper = AmbiguousAlignmentRecordKeeper(
				self.out_prefix, self.gff_dbm, do_overlap_detection=self.do_overlap_detection
			)
		else:
			ambig_bookkeeper = contextlib.nullcontext()

		with ambig_bookkeeper:
			aln_count = 0
			for aln_count, aln in bam_stream:
				start, end, rev_strand = aln.start, aln.end, aln.is_reverse()

				if aln.rid != current_rid:
					# if a new reference sequence is encountered, empty the alignment cache (if needed)
					if current_rid is not None:
						self.process_caches(current_ref)
					current_ref, current_rid = self.bamfile.get_reference(aln.rid)[0], aln.rid

					print("{time}\tNew reference: {ref} ({rid}/{n_ref}). {n_aln} alignments processed.".format(
						time=datetime.now().strftime("%m/%d/%Y,%H:%M:%S"),
						ref=current_ref, rid=aln.rid, n_ref=self.bamfile.n_references(), n_aln=aln_count),
						file=sys.stderr, flush=True)

				if aln.is_ambiguous() and self.require_ambig_bookkeeping():
					# if ambiguous alignments are not treated as individual alignments (dist1 or 1overN mode)
					# then the alignment processing is deferred to the bookkeeper
					ambig_bookkeeper.process_alignment(current_ref, aln, aln_count)
					start, end, rev_strand = None, None, None
				else:
					pair_aligned = all([
						aln.is_paired(),
						aln.rid == aln.rnext,
						not aln.flag & SamFlags.MATE_UNMAPPED
					])
					process_pair = aln.is_unique() or (self.ambig_mode == "primary_only" and aln.is_primary())
					if process_pair and pair_aligned:
						# if a read belongs to a properly-paired (both mates align to the same reference sequence),
						# unique alignment, then it can only be processed if the mate has already been encountered,
						# otherwise it is cached
						start, end, rev_strand = self.umap_cache.process_alignment(aln)

				# at this point only single-end reads, 'improper' and merged pairs should be processed here
				# ( + ambiguous reads in "all1" mode )
				# remember: in all1, pair information is not easy to retain for the secondary alignments!
				if start is not None:
					hits = self.process_alignments_sameref(current_ref, ((current_rid, start, end, rev_strand),))
					if aln.is_unique():
						self.count_manager.update_unique_counts(hits)
					else:
						self.count_manager.update_ambiguous_counts(hits)

			self.process_caches(current_ref)
			t1 = time.time()

			print("Processed {n_align} alignments in {n_seconds:.3f}s.".format(
				n_align=aln_count, n_seconds=t1 - t0), flush=True
			)
			if aln_count == 0:
				print("Warning: bam file does not contain any alignments.")

		if self.require_ambig_bookkeeping():
			return aln_count, ambig_bookkeeper.n_unannotated(), ambig_bookkeeper.dumpfile

		return aln_count, 0, None


	def _process_ambiguous_aln_groups(self, ambig_aln):
		""" iterates over name-sorted (grouped) alignments and processes groups individually """
		print("Processing ambiguous alignments...")
		n_align, current_group = 0, None
		# losing some speed due to the row-iterations here
		# but ok for now
		for aln in ambig_aln.itertuples(index=False, name=None):
			if current_group is None or current_group.qname != aln[0]:
				if current_group:
					n_align += current_group.n_align()
					self.count_manager.update_ambiguous_counts(current_group.resolve())
				current_group = AmbiguousAlignmentGroup(aln)
			else:
				current_group.add_alignment(aln)

		if current_group is not None:
			n_align += current_group.n_align()
			self.count_manager.update_ambiguous_counts(current_group.resolve())

		return n_align


	def process_bamfile(self, bamfile, min_identity=None, min_seqlen=None):
		""" processes one position-sorted bamfile """
		self.bamfile = BamFile(
			bamfile,
			large_header=not self.do_overlap_detection, # ugly!
		)
		# first pass: process uniqs and dump ambigs (if required)
		aln_count, unannotated_ambig, ambig_dumpfile = self.process_alignments(
			min_identity=min_identity,
			min_seqlen=min_seqlen
		)
		self.gff_dbm.clear_caches()

		if aln_count:
			# second pass: process ambiguous alignment groups
			if ambig_dumpfile:
				if os.path.isfile(ambig_dumpfile) and os.stat(ambig_dumpfile).st_size > 0:
					t0 = time.time()

					ambig_aln = FeatureQuantifier.read_ambiguous_alignments(ambig_dumpfile)
					n_align = self._process_ambiguous_aln_groups(ambig_aln)

					t1 = time.time()
					print("Processed {n_align} secondary alignments in {n_seconds:.3f}s.".format(
						n_align=n_align, n_seconds=t1 - t0), flush=True)
				else:
					print(
						"Warning: ambig-mode chosen, "
						"but bamfile does not contain secondary alignments."
					)
					#self.overlap_counter.has_ambig_counts = True  # we expect ambig cols in the outfile!

			self.count_manager.dump_counters(self.out_prefix)
			# to be implemented
			#self.overlap_counter.annotate_counts(
			#	bamfile=self.bamfile,
			#	itermode="counts" if self.do_overlap_detection else "database"
			#)
			#self.overlap_counter.unannotated_reads += unannotated_ambig
			#self.overlap_counter.dump_counts(bam=self.bamfile)

		try:
			if not self.debugmode:
				os.remove(ambig_dumpfile)
		except FileNotFoundError:
			pass

		print("Finished.", flush=True)

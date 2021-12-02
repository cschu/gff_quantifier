import csv
import os
import time
import sys
from datetime import datetime
import contextlib

import pandas

from gffquant.bamreader import BamFile, SamFlags
from gffquant.gff_dbm import GffDatabaseManager
from gffquant.overlap_counter import OverlapCounter
from gffquant.ambig_aln import AmbiguousAlignmentRecordKeeper, AmbiguousAlignmentGroup


class PairedEndAlignmentCache(dict):
	def __init__(self, ambig_alignments=False):
		dict.__init__(self)
		self.ambig_alignments = ambig_alignments

	def empty_cache(self):
		for qname, alignment in self.items():
			n_aln = len(alignment)

			if n_aln == 1:
				# unpaired
				start, end, flag = alignment[0][1:]
				rev_strand = SamFlags.is_reverse_strand(flag)
			elif n_aln == 2:
				# paired
				start, end = BamFile.calculate_fragment_borders(*alignment[0][1:-1], *alignment[1][1:-1])
				rev_strand = None  # TODO: add strand-specific handling by RNAseq protocol
			else:
				raise ValueError(f"{n_aln} primary alignments detected for read {qname}.")

			yield alignment[0][0], start, end, rev_strand

		self.clear()

	def process_alignment(self, alignment):
		# need to keep track of read pairs to avoid dual counts

		start, end, rev_strand = None, None, None

		# check if the mate has already been seen
		mates = self.setdefault(alignment.qname, list())
		if mates:
			if alignment.rnext != mates[0][0]:
				raise ValueError("Alignment ${alignment.qname} seems to be corrupted: {str(alignment)} {str(mates[0])}.")

			# if mate has been seen, calculate the total fragment size and remove the pair from the cache
			start, end = BamFile.calculate_fragment_borders(alignment.start, alignment.end, *mates[0][1:-1])
			del self[alignment.qname]
		else:
			# otherwise cache the first encountered mate and advance to the next read
			mates.append((alignment.rid, alignment.start, alignment.end, alignment.flag))

		return start, end, rev_strand


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
		self.overlap_counter = OverlapCounter(
			out_prefix, self.gff_dbm,
			do_overlap_detection=self.do_overlap_detection, strand_specific=strand_specific,
			feature_distribution=ambig_mode
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
		self.overlap_counter.update_unique_counts(
			self.process_alignments_sameref(ref, self.umap_cache.empty_cache())
		)
		self.overlap_counter.update_ambiguous_counts(
			tuple(self.process_alignments_sameref(ref, self.ambig_cache.empty_cache()))
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
						self.overlap_counter.update_unique_counts(hits)
					else:
						self.overlap_counter.update_ambiguous_counts(tuple(hits))

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
					self.overlap_counter.update_ambiguous_counts(current_group.resolve())
				current_group = AmbiguousAlignmentGroup(aln)
			else:
				current_group.add_alignment(aln)

		if current_group is not None:
			n_align += current_group.n_align()
			self.overlap_counter.update_ambiguous_counts(current_group.resolve())

		print(f"Processed {n_align} secondary alignments.")
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
					self.overlap_counter.has_ambig_counts = True  # we expect ambig cols in the outfile!

			self.overlap_counter.annotate_counts(
				bamfile=self.bamfile,
				itermode="counts" if self.do_overlap_detection else "database"
			)
			self.overlap_counter.unannotated_reads += unannotated_ambig
			self.overlap_counter.dump_counts(bam=self.bamfile)

		try:
			if not self.debugmode:
				os.remove(ambig_dumpfile)
		except FileNotFoundError:
			pass

		print("Finished.", flush=True)

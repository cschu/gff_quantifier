# pylint: disable=C0103

""" module docstring """

import time
import contextlib
import logging

import pandas

from gffquant.bamreader import BamFile, SamFlags
from gffquant.db.annotation_db import AnnotationDatabaseManager
from gffquant.alignment import (
    AmbiguousAlignmentRecordKeeper,
    PairedEndAlignmentCache,
)
from gffquant.counters import CountManager
from gffquant.annotation import DbCountAnnotator, CtCountAnnotator
from gffquant.count_dumper import CountDumper


logger = logging.getLogger(__name__)


class FeatureQuantifier:
    # pylint: disable=R0902,R0913
    TRUE_AMBIG_MODES = ("dist1", "1overN")

    @staticmethod
    def read_ambiguous_alignments(ambig_in):
        """reads the dumped ambig alignments and returns them sorted by alignment group"""
        # using a pandas dataframe/numpy array saves us a lot of memory
        ambig_aln = pandas.read_csv(ambig_in, sep="\t", header=None)
        return ambig_aln.sort_values(axis=0, by=0)

    def __init__(
        self,
        db=None,
        db_index=None,
        out_prefix="gffquant",
        ambig_mode="unique_only",
        reference_type="genome",
        strand_specific=False,
        debugmode=False,
    ):
        self.adm = AnnotationDatabaseManager(db)
        self.umap_cache = PairedEndAlignmentCache()
        self.ambig_cache = PairedEndAlignmentCache(ambig_alignments=True)
        self.do_overlap_detection = reference_type in ("genome", "domain")
        self.count_manager = CountManager(
            distribution_mode=ambig_mode,
            region_counts=reference_type in ("genome", "domain"),
            strand_specific=strand_specific and reference_type not in ("genome", "domain"),
        )
        self.out_prefix = out_prefix
        self.ambig_mode = ambig_mode
        self.bamfile = None
        self.debugmode = debugmode
        self.strand_specific = strand_specific

    def allow_ambiguous_alignments(self):
        """All distribution modes support ambiguous alignments,
        with the exception of `unique_only."""
        return self.ambig_mode != "unique_only"

    def treat_ambiguous_as_unique(self):
        """In 'non-true' ambig modes (all1, primary_only),
        ambiguous alignments are treated as unique alignments."""
        return self.ambig_mode not in FeatureQuantifier.TRUE_AMBIG_MODES

    def require_ambig_bookkeeping(self):
        """'Ambiguous bookkeeping' is a special treatment for true ambig modes.
        It requires all alignments of a read to be processed together."""
        return (
            self.allow_ambiguous_alignments() and not self.treat_ambiguous_as_unique()
        )

    def process_alignments_sameref(self, ref, alignments, aln_count=1):
        """
        Process a group of alignments against the same reference sequence.

        The weird output format is due to compatibility reasons with ambiguous alignments
        when one of the TRUE_AMBIG_MODES is chosen.

        input:
          - reference clear text identifier (required for overlap detection)
          - list of alignments (reference short id, start, end, is_reverse)

        output:
          - generator! with the following items:
            in overlap mode:
              {reference_id:
                  [(begin overlap, end overlap,
				    is_reverse, begin covered region, end covered region)
               for each overlapped region]}
                  + alignment count: 1 if the alignment overlaps any region of interest, else 0
                  + unalignment count: 1 - alignment count
            else:
                  {reference_id:(-1, -1, is_reverse, -1, -1)
                  + alignment count: 1
				  (each reference is already a feature, i.e. everything is aligned)
                  + unalignment count: 0
        """

        for rid, start, end, rev_strand in alignments:
            if self.do_overlap_detection:
                overlaps, coverage = self.adm.get_overlaps(ref, start, end)
                hits = {
                    (ovl.begin, ovl.end, rev_strand, cstart, cend)
                    for ovl, (cstart, cend) in zip(overlaps, coverage)
                }
                # if the alignment overlaps multiple features, each one gets a count
                aln_count = int(bool(overlaps)) * aln_count

            else:
                hits = {(None, None, rev_strand, None, None)}
                # aln_count = 1 / aln_count

            yield ({rid: hits}, aln_count, 0 if aln_count else 1)

    def process_caches(self, ref, aln_count=1):
        """Update the count managers and clear the alignment caches. This is usually done,
        when a new reference is encountered while processing a position-sorted bamfile."""

        self.count_manager.update_counts(
            self.process_alignments_sameref(
                ref, self.umap_cache.empty_cache(), aln_count=aln_count
            ),
            ambiguous_counts=False
        )

        self.count_manager.update_counts(
            self.process_alignments_sameref(
                ref, self.ambig_cache.empty_cache(), aln_count=aln_count
            ),
            ambiguous_counts=True
        )

    def process_alignments(self, ambig_bookkeeper, min_identity=None, min_seqlen=None):
        # pylint: disable=R0914
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
            min_seqlen=min_seqlen,
        )

        aln_count = 0
        for aln_count, aln in bam_stream:
            start, end, rev_strand = aln.start, aln.end, aln.is_reverse()

            if aln.rid != current_rid:
                # if a new reference sequence is encountered,
                # empty the alignment cache (if needed)
                if current_rid is not None:
                    self.process_caches(current_ref)
                current_ref, current_rid = (
                    self.bamfile.get_reference(aln.rid)[0],
                    aln.rid,
                )

            if aln.is_ambiguous() and self.require_ambig_bookkeeping():
                # if ambiguous alignments are not treated as individual alignments
                # (dist1 or 1overN mode)
                # then the alignment processing is deferred to the bookkeeper
                ambig_count = ambig_bookkeeper.get_ambig_alignment_count(aln)
                # ambig_bookkeeper.process_alignment(current_ref, aln, aln_count)
                # start, end, rev_strand = None, None, None
                hits = self.process_alignments_sameref(
                    current_ref, ((current_rid, start, end, rev_strand),), aln_count=ambig_count
                )
                self.count_manager.update_counts(
                    hits, ambiguous_counts=True
                )
                start = None

            else:
                pair_aligned = all(
                    [
                        aln.is_paired(),
                        aln.rid == aln.rnext,
                        not aln.flag & SamFlags.MATE_UNMAPPED,
                    ]
                )
                process_pair = aln.is_unique() or (
                    self.ambig_mode == "primary_only" and aln.is_primary()
                )
                if process_pair and pair_aligned:
                    # if a read belongs to a properly-paired unique alignment
                    # (both mates align to the same reference sequence),
                    # then it can only be processed if the mate
                    # has already been encountered, otherwise it is cached
                    mate = self.umap_cache.process_alignment(aln)
                    mate_start, mate_end, mate_rev_strand = mate
                    if mate_start is not None:
                        hits = self.process_alignments_sameref(
                            current_ref,
                            (
                                (current_rid, start, end, rev_strand),
                                (current_rid, mate_start, mate_end, mate_rev_strand)
                            )
                        )
                        self.count_manager.update_counts(
                            hits, ambiguous_counts=False
                        )

                        mate_start = None

                    # if there was a mate cached, then both mates have been dealt with
                    # hence, signal that the read does not need to be processed further
                    # for overlap/merge later on, the mate_start will be not None
                    start = mate_start

            # at this point only single-end reads, 'improper' and merged pairs
            # should be processed here
            # ( + ambiguous reads in "all1" mode )
            # remember:
            # in all1, pair information is not easy to retain for the secondary alignments!
            if start is not None:
                hits = self.process_alignments_sameref(
                    current_ref, ((current_rid, start, end, rev_strand),)
                )
                self.count_manager.update_counts(
                    hits, ambiguous_counts=not aln.is_unique()
                )

            self.process_caches(current_ref)

        if aln_count == 0:
            print("Warning: bam file does not contain any alignments.")

        t1 = time.time()
        print(f"Processed {aln_count} alignments in {t1 - t0:.3f}s.", flush=True)

        return aln_count, 0, None

    def process_bamfile(self, bamfile, min_identity=None, min_seqlen=None, buffer_size=10000000):
        """processes one position-sorted bamfile"""

        if self.require_ambig_bookkeeping():
            ambig_bookkeeper = AmbiguousAlignmentRecordKeeper(
                self.out_prefix,
                self.adm,
                do_overlap_detection=self.do_overlap_detection,
            )
        else:
            ambig_bookkeeper = contextlib.nullcontext()

        with ambig_bookkeeper:
            self.bamfile = BamFile(
                bamfile,
                large_header=not self.do_overlap_detection,
                buffer_size=buffer_size,
                ambig_bookkeeper=ambig_bookkeeper
            )

            aln_count, unannotated_ambig, _ = self.process_alignments(
                ambig_bookkeeper, min_identity=min_identity, min_seqlen=min_seqlen
            )

            ambig_bookkeeper.clear()

            if aln_count:
                # second pass: process ambiguous alignment groups
                self.count_manager.dump_raw_counters(self.out_prefix, self.bamfile)

                ca_ctr = CtCountAnnotator if self.do_overlap_detection else DbCountAnnotator
                count_annotator = ca_ctr(self.strand_specific)
                count_annotator.annotate(self.bamfile, self.adm, self.count_manager)

                count_dumper = CountDumper(
                    self.out_prefix,
                    has_ambig_counts=self.count_manager.has_ambig_counts(),
                    strand_specific=self.strand_specific,
                )
                count_dumper.dump_feature_counts(
                    self.adm,
                    self.count_manager.get_unannotated_reads() + unannotated_ambig,
                    count_annotator,
                )

                count_dumper.dump_gene_counts(
                    count_annotator.gene_counts,
                    count_annotator.scaling_factors["total_uniq"],
                    count_annotator.scaling_factors["total_ambi"]
                )

        print("Finished.", flush=True)

# pylint: disable=C0103

""" module docstring """

import time
import logging

from gffquant.bamreader import SamFlags
from gffquant.alignment.aln_group import AlignmentGroup
from gffquant.feature_quantifier import FeatureQuantifier
from gffquant.pysam_support import AlignmentProcessor


logger = logging.getLogger(__name__)


class GeneQuantifier(FeatureQuantifier):
    def __init__(
        self,
        db=None,
        out_prefix="gffquant",
        ambig_mode="uniq_only",
        strand_specific=False,
        calc_coverage=False,
    ):
        FeatureQuantifier.__init__(
            self,
            db=db,
            out_prefix=out_prefix,
            ambig_mode=ambig_mode,
            strand_specific=strand_specific,
            reference_type="gene",
            calc_coverage=calc_coverage and False,  # TODO: figure out, but nobody wants it anyway
        )

    def process_bamfile(self, bamfile, min_identity=None, min_seqlen=None, buffer_size=10000000):
        """processes one bamfile"""

        self.alp = AlignmentProcessor(bamfile, "sam")

        # self.bamfile = BamFile(
        #     bamfile,
        #     large_header=not self.do_overlap_detection,
        #     buffer_size=buffer_size
        # )

        aln_count, unannotated_ambig, _ = self.process_alignments(
            min_identity=min_identity, min_seqlen=min_seqlen
        )

        if aln_count:
            self.process_counters(unannotated_ambig)

        print("Finished.", flush=True)

    def process_alignments(self, ambig_bookkeeper=None, min_identity=None, min_seqlen=None):
        # pylint: disable=R0914
        t0 = time.time()

        aln_stream = self.alp.get_alignments(
            min_identity=min_identity,
            min_seqlen=min_seqlen,
            allow_multiple=self.allow_ambiguous_alignments(),
            allow_unique=True,
            filter_flags=SamFlags.SUPPLEMENTARY_ALIGNMENT,
        )

        aln_count = 0
        current_aln_group = None
        for aln_count, aln in enumerate(aln_stream, start=1):
            if self.ambig_mode == "primary_only" and not aln.is_primary():
                continue
            if self.ambig_mode in ("uniq_only", "unique_only") and not aln.is_unique():
                continue

            if current_aln_group is None or current_aln_group.qname != aln.qname:
                if current_aln_group is not None:
                    self.process_alignment_group(current_aln_group)
                current_aln_group = AlignmentGroup()

            current_aln_group.add_alignment(aln)

        if current_aln_group is not None:
            self.process_alignment_group(current_aln_group)

        if aln_count == 0:
            print("Warning: bam file does not contain any alignments.")

        t1 = time.time()
        print(f"Processed {aln_count} alignments in {t1 - t0:.3f}s.", flush=True)

        return aln_count, 0, None

    def process_alignment_group(self, aln_group):
        logging.info("Processing new alignment group %s (%s)", aln_group.qname, aln_group.n_align())
        ambig_counts = aln_group.get_ambig_align_counts()
        if any(ambig_counts) and self.require_ambig_bookkeeping:
            for aln in aln_group.get_alignments():
                if aln is not None:
                    current_ref = self.alp.get_reference(aln.rid)[0]
                    ambig_count = ambig_counts[aln.is_second()]
                    hits = self.process_alignments_sameref(
                        current_ref, (aln.shorten(),), aln_count=ambig_count
                    )
                    self.count_manager.update_counts(
                        hits, ambiguous_counts=True
                    )
        elif aln_group.is_aligned_pair():
            current_ref = self.alp.get_reference(aln_group.primaries[0].rid)[0]
            hits = self.process_alignments_sameref(
                current_ref,
                (
                    aln_group.primaries[0].shorten(),
                    aln_group.primaries[1].shorten(),
                )
            )
            self.count_manager.update_counts(
                hits, ambiguous_counts=False
            )
        else:
            for aln in aln_group.get_alignments():
                current_ref = self.alp.get_reference(aln.rid)[0]
                hits = self.process_alignments_sameref(
                    current_ref, (aln.shorten(),)
                )
                self.count_manager.update_counts(
                    hits, ambiguous_counts=not aln.is_unique()
                )

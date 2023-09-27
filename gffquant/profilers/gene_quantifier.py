# pylint: disable=C0103

""" module docstring """

import logging

from .feature_quantifier import FeatureQuantifier
from .. import __tool__

logger = logging.getLogger(__name__)


class GeneQuantifier(FeatureQuantifier):
    # pylint: disable=R0913,W0511,R0801
    def __init__(
        self,
        db=None,
        out_prefix=__tool__,
        ambig_mode="uniq_only",
        strand_specific=False,
        paired_end_count=1,
    ):
        FeatureQuantifier.__init__(
            self,
            db=db,
            out_prefix=out_prefix,
            ambig_mode=ambig_mode,
            strand_specific=strand_specific,
            reference_type="gene",
            paired_end_count=paired_end_count,
        )

    def process_alignment_group(self, aln_group, aln_reader):
        # logger.info("Processing new alignment group %s (%s)", aln_group.qname, aln_group.n_align())
        ambig_counts = aln_group.get_ambig_align_counts()
        if any(ambig_counts) and self.require_ambig_bookkeeping:
            for aln in aln_group.get_alignments():
                current_ref = self.register_reference(aln.rid, aln_reader)
                ambig_count = ambig_counts[aln.is_second()]
                hits = self.process_alignments_sameref(
                    current_ref, (aln.shorten(),), aln_count=ambig_count
                )
                self.count_manager.update_counts(
                    hits, ambiguous_counts=True, pair=aln_group.is_paired(), pe_library=aln_group.pe_library,
                )
        elif aln_group.is_aligned_pair():
            current_ref = self.register_reference(aln_group.primaries[0].rid, aln_reader)
            hits = self.process_alignments_sameref(
                current_ref,
                (
                    aln_group.primaries[0].shorten(),
                    aln_group.primaries[1].shorten(),
                )
            )
            self.count_manager.update_counts(
                hits, ambiguous_counts=False, pair=True, pe_library=aln_group.pe_library,
            )
        else:
            for aln in aln_group.get_alignments():
                current_ref = self.register_reference(aln.rid, aln_reader)
                hits = self.process_alignments_sameref(
                    current_ref, (aln.shorten(),)
                )
                self.count_manager.update_counts(
                    hits,
                    ambiguous_counts=not aln.is_unique(),
                    pair=aln_group.is_paired(),
                    pe_library=aln_group.pe_library,
                )

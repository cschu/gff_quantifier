# pylint: disable=C0103

""" module docstring """

from itertools import count
import logging

from gffquant.db.annotation_db import AnnotationDatabaseManager
from gffquant.counters import CountManager
from gffquant.annotation import GeneCountAnnotator, RegionCountAnnotator
from gffquant.count_dumper import CountDumper
from gffquant.coverage_counter import CoverageCounter


logger = logging.getLogger(__name__)


class FeatureQuantifier:
    # pylint: disable=R0902,R0913
    TRUE_AMBIG_MODES = ("dist1", "1overN")

    def __init__(
        self,
        db=None,
        out_prefix="gffquant",
        ambig_mode="unique_only",
        reference_type="genome",
        strand_specific=False,
        calc_coverage=False,
    ):
        self.db = db
        self.adm = None
        self.do_overlap_detection = reference_type in ("genome", "domain")
        self.count_manager = CountManager(
            distribution_mode=ambig_mode,
            region_counts=reference_type in ("genome", "domain"),
            strand_specific=strand_specific and reference_type not in ("genome", "domain"),
            calc_coverage=calc_coverage,
        )
        self.out_prefix = out_prefix
        self.ambig_mode = ambig_mode
        self.bamfile = None
        self.alp = None
        self.strand_specific = strand_specific
        self.calc_coverage = calc_coverage

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

    def process_counters(self, unannotated_ambig):
        if self.adm is None:
            self.adm = AnnotationDatabaseManager(self.db)

        self.count_manager.dump_raw_counters(self.out_prefix, self.alp)

        cov_ctr = CoverageCounter() if self.calc_coverage else None
        
        if self.do_overlap_detection:
            count_annotator = RegionCountAnnotator(self.strand_specific)
            count_annotator.annotate(self.alp, self.adm, self.count_manager, coverage_counter=cov_ctr)
        else:
            count_annotator = GeneCountAnnotator(self.strand_specific)
            count_annotator.annotate(self.alp, self.adm, self.count_manager)

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

        if self.calc_coverage:
            print(*cov_ctr.items(), sep="\n")
            count_dumper.dump_coverage(self.adm, cov_ctr)

        self.adm.clear_caches()

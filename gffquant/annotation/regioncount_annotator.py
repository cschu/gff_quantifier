import numpy as np

from . import CountAnnotator
from ..counters import AlignmentCounter


class RegionCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for contig/region-based counting. """

    def __init__(self, strand_specific, report_scaling_factors=True):
        CountAnnotator.__init__(self, strand_specific, report_scaling_factors=report_scaling_factors)

    # pylint: disable=R0914,W0613
    def annotate(self, refmgr, db, counter: AlignmentCounter, gene_group_db=False):
        """
        Annotate a set of region counts via db-lookup.
        input:
        - bam: bamr.BamFile to use as lookup table for reference names
        - db: GffDatabaseManager holding functional annotation database
        """
        for rid in counter.get_all_regions(region_counts=True):
            ref = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)[0]

            for region in counter.get_regions(rid):
                if self.strand_specific:
                    (start, end), rev_strand = region
                else:
                    (start, end), rev_strand = region, None
                # the region_annotation is a tuple of key-value pairs:
                # (strand, func_category1: subcategories, func_category2: subcategories, ...)
                # the first is the strand, the second is the gene id, the rest are the features

                region_annotation = db.query_sequence(ref, start=start, end=end)
                if region_annotation is not None:
                    region_strand, feature_id, region_annotation = region_annotation
                    if feature_id is None:
                        feature_id = ref

                    on_other_strand = (region_strand == "+" and rev_strand) \
                        or (region_strand == "-" and not rev_strand)

                    antisense_region = self.strand_specific and on_other_strand

                    uniq_counts, ambig_counts = counter.get_counts(
                        (rid, start, end), region_counts=True, strand_specific=self.strand_specific
                    )

                    if self.strand_specific:
                        # if the region is antisense, 'sense-counts' (relative to the) region come from the
                        # negative strand and 'antisense-counts' from the positive strand
                        # vice-versa for a sense-region
                        strand_specific_counts = (
                            (counter.MINUS_STRAND, counter.PLUS_STRAND)
                            if antisense_region
                            else (counter.PLUS_STRAND, counter.MINUS_STRAND)
                        )
                    else:
                        strand_specific_counts = None

                    region_length = end - start + 1
                    counts = self.compute_count_vector(
                        uniq_counts,
                        ambig_counts,
                        region_length,
                        strand_specific_counts=strand_specific_counts,
                        region_counts=True,
                    )

                    self.distribute_feature_counts(counts, region_annotation)

                    gcounts = self.gene_counts.setdefault(
                        feature_id, np.zeros(self.bins)
                    )
                    gcounts += counts
                    self.total_gene_counts += counts[:4]

        self.calculate_scaling_factors()

# pylint: disable=C0103,W0105,R0902

""" This module contains code for transforming gene counts to feature counts. """

import logging

from itertools import chain

import numpy as np


logger = logging.getLogger(__name__)

"""
normalizeCounts nmethod counts sizes
	| nmethod `elem` [NMScaled, NMFpkm] = do
		-- count vectors always include a -1 at this point (it is
		-- ignored in output if the user does not request it, but is
		-- always computed). Thus, we compute the sum without it and do
		-- not normalize it later:
		let totalCounts v = withVector v (VU.sum . VU.tail)
		initial <- totalCounts counts
		normalizeCounts NMNormed counts sizes
		afternorm <- totalCounts counts
		let factor
				| nmethod == NMScaled = initial / afternorm
				| otherwise = 1.0e9 / initial --- 1e6 [million fragments] * 1e3 [kilo basepairs] = 1e9
		liftIO $ forM_ [1.. VUM.length counts - 1] (VUM.unsafeModify counts (* factor))
"""


class CountAnnotator(dict):
    """ CountAnnotator is the parent class for the two different count annotators. """

    def __init__(self, strand_specific, report_scaling_factors=True):
        """
        input:
        - strand_specific: true | false
        """
        dict.__init__(self)
        self.strand_specific = strand_specific
        self.bins = 12 if strand_specific else 4

        # [uniq_raw, uniq_normed, combined_raw, combined_normed]
        self.total_counts = np.zeros(4)
        self.total_gene_counts = np.zeros(4)
        self.unannotated_counts = np.zeros(4)
        # holds total_counts-like vectors for feature-wise scaling factor calculation
        self.feature_count_sums = {}
        self.scaling_factors = {}
        self.gene_counts = {}
        self.report_scaling_factors = report_scaling_factors

    def distribute_feature_counts(self, counts, region_annotation):
        """
        Distributes the counts for a region/gene among the annotated functional
        categories/features

        input:
        - counts: [uniq_raw, uniq_normed, combined_raw, combined_normed]
        - region_annotation: functional categories/features
        """

        # add region-counts once for gene/region scaling factor calculation
        self.total_counts += counts[:4]

        for category, category_counts in region_annotation:
            # add category-counts once for category scaling factor calculation
            total_fcounts = self.feature_count_sums.setdefault(category, np.zeros(4))
            total_fcounts += counts[:4]

            for feature in category_counts:
                self.add_counts(category, feature, counts)

            if category_counts:
                # eggnog-mapper annotation tables in v2.1.2 (maybe earlier?) have '-' instead of empty cells
                # category counts may be empty in case a non-patched db based on such tables is used
                # in that case we're adding counts of non-existing features to the category
                # without checking for empty cat counts!
                self.add_counts(category, f"cat:::{category}", counts)

    def add_counts(self, category, feature, counts):
        """ Increments feature counts by input count vector """
        fcounts = self.setdefault(category, {}).setdefault(feature, np.zeros(self.bins))
        fcounts += counts

    def calculate_scaling_factors(self, default_scaling_factor=0):
        """ Calculates all scaling factors.
        scaling_factor = uniq_counts / normed_counts
        input:
        - default_scaling_factor: if normed counts do not exist, return this instead
        """

        def calc_scaling_factor(raw, normed, default=0):
            return (raw / normed) if normed else default

        total_uniq, total_uniq_normed, total_ambi, total_ambi_normed = self.total_counts
        logger.info(
            "TOTAL COUNTS: uraw=%s unorm=%s araw=%s anorm=%s",
            total_uniq, total_uniq_normed, total_ambi, total_ambi_normed
        )

        self.scaling_factors["total_uniq"] = calc_scaling_factor(
            total_uniq, total_uniq_normed, default_scaling_factor
        )

        self.scaling_factors["total_ambi"] = calc_scaling_factor(
            total_ambi, total_ambi_normed, default_scaling_factor
        )

        total_uniq, total_uniq_normed, total_ambi, total_ambi_normed = self.total_gene_counts
        logger.info(
            "TOTAL GENE COUNTS: uraw=%s unorm=%s araw=%s anorm=%s",
            total_uniq, total_uniq_normed, total_ambi, total_ambi_normed
        )

        self.scaling_factors["total_gene_uniq"] = calc_scaling_factor(
            total_uniq, total_uniq_normed, default_scaling_factor
        )

        self.scaling_factors["total_gene_ambi"] = calc_scaling_factor(
            total_ambi, total_ambi_normed, default_scaling_factor
        )

        fc_items = self.feature_count_sums.items()
        for category, (
            total_uniq,
            total_uniq_normed,
            total_ambi,
            total_ambi_normed,
        ) in fc_items:

            self.scaling_factors[category] = (
                calc_scaling_factor(
                    total_uniq, total_uniq_normed, default_scaling_factor
                ),
                calc_scaling_factor(
                    total_ambi, total_ambi_normed, default_scaling_factor
                )
            )

            if self.report_scaling_factors:
                logger.info(
                    "Calculating scaling factors for category=%s: uraw=%s unorm=%s araw=%s anorm=%s -> factors=%s",
                    category, total_uniq, total_uniq_normed,
                    total_ambi, total_ambi_normed, self.scaling_factors[category]
                )

    # pylint: disable=R0913
    def compute_count_vector(
        self,
        uniq_counts,
        ambig_counts,
        length,
        strand_specific_counts=None,
        region_counts=False,
    ):
        """Computes a count vector for a region."""
        # we have either 4 bins (unstranded) or 12 (strand-specific)
        # UNSTRANDED = {uniq,ambig} x {raw,normalised}
        # STRANDED = UNSTRANDED x {all,sense/plus,antisense/minus}
        counts = np.zeros(self.bins)

        if strand_specific_counts is not None:
            ss_counts, as_counts = strand_specific_counts
            # counts[4:12] are strand-specific values
            # uniq_raw, uniq_norm, combined_raw, combined_norm
            counts[4:8] = uniq_counts[ss_counts]
            counts[8:12] = uniq_counts[as_counts]
            # add the ambig counts to combined_raw, combined_norm
            counts[6:8] += ambig_counts[ss_counts]
            counts[10:12] += ambig_counts[as_counts]

        # the first 4 elements (counts[0:4]) are unstranded:
        # uniq_raw, uniq_norm, combined_raw, combined_norm
        # 1. each of these fields gets a copy of the unique count sum
        # 2. add the ambiguous counts to the combined_ elements

        # pylint: disable=R1727
        if region_counts and False:
            # used to ask for region_counts and coverage_counts
            counts[0:4] = sum(x[2] for x in chain(*uniq_counts) if x is not None)
            counts[2:4] += sum(x[2] for x in chain(*ambig_counts) if x is not None)
        else:
            counts[0:4] = sum(uniq_counts)
            counts[2:4] += sum(ambig_counts)

        # 3. all odd elements (including strand-specific) are length-normalised
        counts[1::2] /= float(length)

        return counts


class RegionCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for contig/region-based counting. """

    def __init__(self, strand_specific, report_scaling_factors=True):
        CountAnnotator.__init__(self, strand_specific, report_scaling_factors=report_scaling_factors)

    # pylint: disable=R0914
    def annotate(self, refmgr, db, count_manager):
        """
        Annotate a set of region counts via db-lookup.
        input:
        - bam: bamr.BamFile to use as lookup table for reference names
        - db: GffDatabaseManager holding functional annotation database
        - count_manager: count_data
        """
        for rid in set(count_manager.uniq_regioncounts).union(
            count_manager.ambig_regioncounts
        ):
            ref = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)[0]

            for region in count_manager.get_regions(rid):
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

                    uniq_counts, ambig_counts = count_manager.get_counts(
                        (rid, start, end), region_counts=True, strand_specific=self.strand_specific
                    )

                    if self.strand_specific:
                        # if the region is antisense, 'sense-counts' (relative to the) region come from the
                        # negative strand and 'antisense-counts' from the positive strand
                        # vice-versa for a sense-region
                        strand_specific_counts = (
                            (count_manager.MINUS_STRAND, count_manager.PLUS_STRAND)
                            if antisense_region
                            else (count_manager.PLUS_STRAND, count_manager.MINUS_STRAND)
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


class GeneCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for gene-based counting. """

    def __init__(self, strand_specific, report_scaling_factors=True):
        CountAnnotator.__init__(self, strand_specific, report_scaling_factors=report_scaling_factors)

    def annotate(self, refmgr, db, count_manager):
        """
        Annotate a set of gene counts via db-iteration.
        input:
        - bam: bamr.BamFile to use as reverse lookup table for reference ids
        - db: GffDatabaseManager holding functional annotation database
        - count_manager: count_data
        """
        strand_specific_counts = (
            (count_manager.PLUS_STRAND, count_manager.MINUS_STRAND)
            if self.strand_specific else None
        )

        for rid in set(count_manager.uniq_seqcounts).union(
            count_manager.ambig_seqcounts
        ):
            ref, region_length = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)

            uniq_counts, ambig_counts = count_manager.get_counts(
                rid, region_counts=False, strand_specific=self.strand_specific
            )

            counts = self.compute_count_vector(
                uniq_counts,
                ambig_counts,
                region_length,
                strand_specific_counts=strand_specific_counts,
            )

            gcounts = self.gene_counts.setdefault(ref, np.zeros(self.bins))
            gcounts += counts
            self.total_gene_counts += counts[:4]

            region_annotation = db.query_sequence(ref)
            if region_annotation is not None:
                _, _, region_annotation = region_annotation
                self.distribute_feature_counts(counts, region_annotation)
            else:
                self.unannotated_counts += counts[:4]

        self.calculate_scaling_factors()

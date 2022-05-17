# pylint: disable=C0103

import numpy as np


""" This module contains code for transforming gene counts to feature counts. """


class CountAnnotator(dict):
    """ CountAnnotator is the parent class for the two different count annotators. """

    def __init__(self, strand_specific):
        """
        input:
        - strand_specific: true | false
        """
        dict.__init__(self)
        self.strand_specific = strand_specific
        self.bins = 12 if strand_specific else 4

        # [uniq_raw, uniq_normed, combined_raw, combined_normed]
        self.total_counts = np.zeros(4)
        # holds total_counts-like vectors for feature-wise scaling factor calculation
        self.feature_count_sums = {}
        self.scaling_factors = {}
        self.gene_counts = {}

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
            print("DFC:", category, self.total_counts, total_fcounts)

            for feature in category_counts:
                self.add_counts(category, feature, counts)
                print("DFC:", category, feature, counts)

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

        print("TOTAL COUNTS:", self.total_counts)
        total_uniq, total_uniq_normed, total_ambi, total_ambi_normed = self.total_counts

        self.scaling_factors["total_uniq"] = calc_scaling_factor(
            total_uniq, total_uniq_normed, default_scaling_factor
        )

        self.scaling_factors["total_ambi"] = calc_scaling_factor(
            total_ambi, total_uniq_normed, default_scaling_factor
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

            print(
                "SCALING FACTOR", category, total_uniq, total_uniq_normed,
                total_ambi, total_ambi_normed, self.scaling_factors[category]
            )

    def compute_count_vector(self, count_manager, rid, length, antisense_region=False):
        """Computes a count vector for a region."""
        # we have either 4 bins (unstranded) or 12 (strand-specific)
        # UNSTRANDED = {uniq,ambig} x {raw,normalised}
        # STRANDED = UNSTRANDED x {all,sense/plus,antisense/minus}
        counts = np.zeros(self.bins)

        # not sure if this is still needed
        # use_region_counts = isinstance(self, CtCountAnnotator)

        uniq_counts, ambig_counts = count_manager.get_counts(
            rid, region_counts=False, strand_specific=self.strand_specific
        )

        if self.strand_specific:
            # if the region is antisense, 'sense-counts' (relative to the) region come from the
            # negative strand and 'antisense-counts' from the positive strand
            # vice-versa for a sense-region
            ss_counts, as_counts = (
                (count_manager.MINUS_STRAND, count_manager.PLUS_STRAND)
                if antisense_region
                else (count_manager.PLUS_STRAND, count_manager.MINUS_STRAND)
            )

            # counts[4:12] are strand-specific values
            # uniq_raw, uniq_norm, combined_raw, combined_norm
            counts[4:8] = uniq_counts[ss_counts]
            counts[8:12] = uniq_counts[as_counts]
            # add the ambig counts to combined_raw, combined_norm
            counts[6:8] += ambig_counts[ss_counts]
            counts[10:12] += ambig_counts[as_counts]

        # counts[0:4] are unstranded
        # uniq_raw, uniq_norm, combined_raw, combined_norm
        counts[0:4] = sum(uniq_counts)
        # add the ambig counts to combined_raw, combined_norm
        counts[2:4] += sum(ambig_counts)
        # all odd elements are length-normalised
        counts[1::2] /= float(length)

        return counts


class CtCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for contig/region-based counting. """

    def __init__(self, strand_specific):
        CountAnnotator.__init__(self, strand_specific)

    def annotate(self, bam, db, count_manager):
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
            ref = bam.get_reference(rid[0] if isinstance(rid, tuple) else rid)[0]
            for start, end, rev_strand in count_manager.get_regions(rid):
                # region_annotation is a tuple of key-value pairs:
                # (strand, func_category1: subcategories, func_category2: subcategories, ...)
                # the first is the strand, the second is the gene id, the rest are the features
                region_annotation = db.get_data(ref, start, end)
                region_strand = region_annotation[0][1]

                on_other_strand = (region_strand == "+" and rev_strand) \
                    or (region_strand == "-" and not rev_strand)

                antisense_region = self.strand_specific and on_other_strand

                counts = self.compute_count_vector(
                    count_manager,
                    (rid, start, end),
                    end - start + 1,
                    antisense_region=antisense_region,
                )

                # this is [2:] as the gene id comes in the feature string as described
                # above! (region key != gene id!!)
                self.distribute_feature_counts(counts, region_annotation[2:])

                gcounts = self.gene_counts.setdefault(
                    region_annotation[1][1], np.zeros(self.bins)
                )
                gcounts += counts

        self.calculate_scaling_factors()


class DbCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for gene-based counting. """

    def __init__(self, strand_specific):
        CountAnnotator.__init__(self, strand_specific)

    def annotate(self, bam, db, count_manager):
        """
        Annotate a set of gene counts via db-iteration.
        input:
        - bam: bamr.BamFile to use as reverse lookup table for reference ids
        - db: GffDatabaseManager holding functional annotation database
        - count_manager: count_data
        """
        for ref, region_annotation in db.iterate():
            rid = bam.revlookup_reference(ref)
            if rid is not None:
                _, region_length = bam.get_reference(rid[0] if isinstance(rid, tuple) else rid)
                counts = self.compute_count_vector(count_manager, rid, region_length)
                self.distribute_feature_counts(counts, region_annotation[1:])

                gcounts = self.gene_counts.setdefault(ref, np.zeros(self.bins))
                gcounts += counts

        self.calculate_scaling_factors()

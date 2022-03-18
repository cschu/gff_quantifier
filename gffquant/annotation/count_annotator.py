# pylint: disable=C0103

""" module docstring """

import numpy as np


class CountAnnotator(dict):
    def __init__(self, strand_specific):
        dict.__init__(self)
        self.strand_specific = strand_specific
        self.bins = 12 if strand_specific else 4

        self.total_counts = np.zeros(4)
        self.feature_count_sums = {}
        self.scaling_factors = {}
        self.gene_counts = {}

    def distribute_feature_counts(self, counts, region_annotation):
        for category, category_counts in region_annotation:
            total_fcounts = self.feature_count_sums.setdefault(category, np.zeros(4))

            for feature in category_counts:
                self.add_counts(category, feature, counts)
                print("DFC:", category, feature, counts)

            self.total_counts += counts[:4]
            total_fcounts += counts[:4]
            print("DFC:", category, self.total_counts, total_fcounts)

    def add_counts(self, category, feature, counts):
        fcounts = self.setdefault(category, {}).setdefault(feature, np.zeros(self.bins))
        fcounts += counts

    def calculate_scaling_factors(self, default_scaling_factor=0):
        total_uniq, total_uniq_normed, total_ambi, total_ambi_normed = self.total_counts

        self.scaling_factors["total_uniq"] = (
            (total_uniq / total_uniq_normed)
            if total_uniq_normed
            else default_scaling_factor
        )
        self.scaling_factors["total_ambi"] = (
            (total_ambi / total_ambi_normed)
            if total_ambi_normed
            else default_scaling_factor
        )

        fc_items = self.feature_count_sums.items()
        for category, (
            total_uniq,
            total_uniq_normed,
            total_ambi,
            total_ambi_normed,
        ) in fc_items:

            self.scaling_factors[category] = (
                (total_uniq / total_uniq_normed)
                if total_uniq_normed
                else default_scaling_factor,
                (total_ambi / total_ambi_normed)
                if total_ambi_normed
                else default_scaling_factor,
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
    def __init__(self, strand_specific):
        CountAnnotator.__init__(self, strand_specific)

    def annotate(self, bam, db, count_manager):
        for rid in set(count_manager.uniq_regioncounts).union(
            count_manager.ambig_regioncounts
        ):
            ref = bam.get_reference(rid)[0]
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

                self.distribute_feature_counts(counts, region_annotation[2:])

                gcounts = self.gene_counts.setdefault(
                    region_annotation[1][1], np.zeros(self.bins)
                )
                gcounts += counts

        self.calculate_scaling_factors()


class DbCountAnnotator(CountAnnotator):
    def __init__(self, strand_specific):
        CountAnnotator.__init__(self, strand_specific)

    def annotate(self, bam, db, count_manager):

        for ref, region_annotation in db.iterate():
            rid = bam.revlookup_reference(ref)
            if rid is not None:
                _, region_length = bam.get_reference(rid)
                counts = self.compute_count_vector(count_manager, rid, region_length)
                self.distribute_feature_counts(counts, region_annotation[1:])

                gcounts = self.gene_counts.setdefault(ref, np.zeros(self.bins))
                gcounts += counts

        self.calculate_scaling_factors()

# pylint: disable=R0914

""" module docstring """
import logging

import numpy as np

from .count_annotator import CountAnnotator
from ..counters import AlignmentCounter
from ..counters.count_matrix import CountMatrix
from ..db.annotation_db import AnnotationDatabaseManager


logger = logging.getLogger(__name__)


class GeneCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for gene-based counting. """
    def __init__(self, strand_specific, report_scaling_factors=True):
        """ __init__() """
        CountAnnotator.__init__(self, strand_specific, report_scaling_factors=report_scaling_factors)

    def annotate_gene_counts(
        self,
        refmgr,
        db: AnnotationDatabaseManager,
        counter: AlignmentCounter,
        filtered_reads: float,
        gene_group_db=False
    ):
        categories = list(db.get_categories())
        category_sums = np.zeros((len(categories), 8))
        functional_counts = CountMatrix(8)
        unannotated_counts = CountMatrix(8, 1)

        for rid, counts in counter:
            if gene_group_db:
                if rid == "0":
                    unannotated_counts[CountMatrix.NO_ANNOTATION] += counts
                    continue
                ggroup_id = int(rid, 16)
            else:
                ref, _ = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)
                ggroup_id = ref

            region_annotation = db.query_sequence(ggroup_id, grouped_db=gene_group_db,)
            if region_annotation is not None:
                *_, region_annotation = region_annotation
                for category_id, features in region_annotation:
                    category_id = int(category_id)
                    category_sums[category_id] += counts
                    for feature_id in features:
                        feature_id = int(feature_id)
                        functional_counts[(category_id, feature_id)] += counts
            elif not gene_group_db:
                unannotated_counts[CountMatrix.NO_ANNOTATION] += counts

        functional_counts.drop_unindexed()

        for i, category in enumerate(categories):
            u_sf, c_sf = (
                CountMatrix.calculate_scaling_factor(*category_sums[i][CountMatrix.RAW_COLUMNS[0]:CountMatrix.LNORM_COLUMNS[0] + 1:]),
                CountMatrix.calculate_scaling_factor(*category_sums[i][CountMatrix.RAW_COLUMNS[1]:CountMatrix.LNORM_COLUMNS[1] + 1:]),
            )
            rpkm_sf = 1e9 / filtered_reads

            rows = tuple(
                key[0] == category.id
                for key, _ in functional_counts
            )

            functional_counts.scale_column(CountMatrix.LNORM_COLUMNS[0], CountMatrix.SCALED_COLUMNS[0], u_sf, rows=rows)
            functional_counts.scale_column(CountMatrix.LNORM_COLUMNS[1], CountMatrix.SCALED_COLUMNS[1], c_sf, rows=rows)
            functional_counts.scale_column(CountMatrix.LNORM_COLUMNS[0], CountMatrix.SCALED_COLUMNS[0], rpkm_sf, rows=rows)
            functional_counts.scale_column(CountMatrix.LNORM_COLUMNS[1], CountMatrix.SCALED_COLUMNS[1], rpkm_sf, rows=rows)

            category_sums[i, CountMatrix.SCALED_COLUMNS[0]] = category_sums[i, CountMatrix.LNORM_COLUMNS[0]] * u_sf
            category_sums[i, CountMatrix.SCALED_COLUMNS[1]] = category_sums[i, CountMatrix.LNORM_COLUMNS[1]] * c_sf
            category_sums[i, CountMatrix.RPKM_COLUMNS[0]] = category_sums[i, CountMatrix.LNORM_COLUMNS[0]] * rpkm_sf
            category_sums[i, CountMatrix.RPKM_COLUMNS[1]] = category_sums[i, CountMatrix.LNORM_COLUMNS[1]] * rpkm_sf

        return functional_counts, category_sums, unannotated_counts

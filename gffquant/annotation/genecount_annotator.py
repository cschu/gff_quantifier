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
        gene_group_db=False
    ):
        categories = list(db.get_categories())
        category_sums = np.zeros((len(categories), 6))
        functional_counts = CountMatrix(6)

        # for category in categories:
        #     features = ((feature.name, feature) for feature in db.get_features(category.id))
        #     for _, feature in sorted(features, key=lambda x: x[0]):
        #         _ = functional_counts[(category.id, feature.id)]

        for rid, counts in counter:
            counts = counter[rid]
            if gene_group_db:
                ggroup_id = rid
            else:
                ref, _ = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)
                ggroup_id = ref

            region_annotation = db.query_sequence(ggroup_id)
            if region_annotation is not None:
                _, _, region_annotation = region_annotation
                for category_id, features in region_annotation:
                    category_id = int(category_id)
                    category_sums[category_id] += counts
                    for feature_id in features:
                        feature_id = int(feature_id)
                        functional_counts[(category_id, feature_id)] += counts

        functional_counts.drop_unindexed()

        for i, category in enumerate(categories):
            u_sf, c_sf = (
                CountMatrix.calculate_scaling_factor(*category_sums[i][0:2]),
                CountMatrix.calculate_scaling_factor(*category_sums[i][3:5]),
            )

            rows = tuple(
                key[0] == category.id
                for key, _ in functional_counts
            )

            functional_counts.scale_column(1, u_sf, rows=rows)
            functional_counts.scale_column(4, c_sf, rows=rows)

            category_sums[i, 2] = category_sums[i, 1] * u_sf
            category_sums[i, 5] = category_sums[i, 4] * c_sf

        return functional_counts, category_sums

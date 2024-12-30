# pylint: disable=R0914

""" module docstring """
import logging

import numpy as np

from .count_annotator import CountAnnotator
from ..counters import AlignmentCounter
from ..db.annotation_db import AnnotationDatabaseManager


logger = logging.getLogger(__name__)


class GeneCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for gene-based counting. """

    def __init__(self, strand_specific, report_scaling_factors=True):
        """ __init__() """
        CountAnnotator.__init__(self, strand_specific, report_scaling_factors=report_scaling_factors)

    def annotate2(self, refmgr, db: AnnotationDatabaseManager, counter: AlignmentCounter, gene_group_db=False):
        for category in db.get_categories():
            features = tuple(db.get_features(category.id))
            category_counts = np.zeros(
                (len(features) + 1, 6,),
                dtype='float64',
            )
            category_index = {
                feature.id: i
                for i, feature in enumerate(features, start=1)
            }
            category_names = {
                feature.id: feature.name
                for feature in features
            }
            for rid in counter:
                # counts = counter.get_counts(rid, strand_specific=self.strand_specific)
                counts = counter[rid]
                if gene_group_db:
                    # gene_id, ggroup_id = rid, rid
                    ggroup_id = rid
                else:
                    ref, _ = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)
                    # gene_id, ggroup_id = ref, ref
                    ggroup_id = ref

                region_annotation = db.query_sequence(ggroup_id)
                # logger.info("REGION_ANNOTATION: %s (%s)", str(region_annotation), ggroup_id)
                if region_annotation is not None:
                    _, _, region_annotation = region_annotation
                    category_features = dict(region_annotation).get(str(category.id))
                    if category_features is not None:
                        category_counts[0] += counts  # category row
                        for cf in category_features:
                            category_counts[category_index.get(int(cf))] += counts

                # elif it == 0:
                #     self.unannotated_counts += counts[:4]

            # count_sums = category_counts[1:].sum(axis=0)
            count_sums = category_counts[0]

            # should scaled counts use a factor derived from all counts
            # or should multi-feature counts only contribute once?
            # uniq_scaling_factor = (count_sums[0] / count_sums[1], 1.0)[count_sums[1] == 0]
            # ambig_scaling_factor = (count_sums[2] / count_sums[3], 1.0)[count_sums[3] == 0]
            # pre 2.19 category count scaling was based on total counts
            # uniq_scaling_factor, ambig_scaling_factor = 1.0, 1.0
            # if category_counts[0][1]:
            #     uniq_scaling_factor = category_counts[0][0] / category_counts[0][1]
            # if category_counts[0][3]:
            #     ambig_scaling_factor = category_counts[0][2] / category_counts[0][3]
            uniq_scaling_factor, combined_scaling_factor = (
                AlignmentCounter.calculate_scaling_factor(*count_sums[0:2]),
                AlignmentCounter.calculate_scaling_factor(*count_sums[3:5]),
            )

            # apply scaling factors
            category_counts[:, 2] = category_counts[:, 1] * uniq_scaling_factor
            category_counts[:, 5] = category_counts[:, 4] * combined_scaling_factor

            logger.info(
                "GCA:: %s CATEGORY COUNTS: uraw=%s unorm=%s araw=%s anorm=%s => SF: %s %s",
                category.name,
                # *count_sums,
                count_sums[0], count_sums[1], count_sums[3], count_sums[4],
                uniq_scaling_factor, combined_scaling_factor,
            )

            yield (
                category.name,
                category_counts,
                category_index,
                category_names,
                uniq_scaling_factor,
                combined_scaling_factor,
            )

    def annotate(self, refmgr, db: AnnotationDatabaseManager, counter: AlignmentCounter, gene_group_db=False):
        """ Annotate a set of gene counts with functional annotations. """

        # formerly used in compute_count_vector
        strand_specific_counts = (
            (counter.PLUS_STRAND, counter.MINUS_STRAND)
            if self.strand_specific else None
        )

        for rid in counter.get_all_regions():
            counts = counter.get_counts(rid, strand_specific=self.strand_specific)

            if gene_group_db:
                # ref_tokens = ref.split(".")
                # gene_id, ggroup_id = ".".join(ref_tokens[:-1]), ref_tokens[-1]
                # gene_id, ggroup_id = rid, rid
                ggroup_id = rid
            else:
                ref, _ = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)
                # gene_id, ggroup_id = ref, ref
                ggroup_id = ref

            region_annotation = db.query_sequence(ggroup_id)
            if region_annotation is not None:
                _, _, region_annotation = region_annotation
                self.distribute_feature_counts(counts, region_annotation)

            else:
                # logger.info("GCAnnotator: Gene %s (group=%s) has no information in database.", gene_id, ggroup_id)
                self.unannotated_counts += counts[:4]

        self.calculate_scaling_factors()

""" module docstring """
import logging

import numpy as np

from .count_annotator import CountAnnotator
from .count_writer import CountWriter
from ..counters import AlignmentCounter
from ..db.annotation_db import AnnotationDatabaseManager


logger = logging.getLogger(__name__)


class GeneCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for gene-based counting. """

    def __init__(self, strand_specific, report_scaling_factors=True):
        """ __init__() """
        CountAnnotator.__init__(self, strand_specific, report_scaling_factors=report_scaling_factors)

    def annotate2(self, refmgr, db: AnnotationDatabaseManager, counter: AlignmentCounter, gene_group_db=False):
        for it, category in enumerate(db.get_categories()):
            features = tuple(db.get_features(category.id))
            # total_reads     483808.00000    483808.00000    483808.00000    483808.00000    483808.00000    483808.00000
            # filtered_reads  454437.00000    454437.00000    454437.00000    454437.00000    454437.00000    454437.00000
            # category        45359.50000     47.10706        42266.81963     152875.83896    224.72779       149853.25971
            category_counts = np.zeros(
                (len(features) + 1, 4,),
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
            for rid in counter.get_all_regions():
                counts = counter.get_counts(rid, strand_specific=self.strand_specific)
                if gene_group_db:
                    gene_id, ggroup_id = rid, rid
                else:
                    ref, _ = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)
                    gene_id, ggroup_id = ref, ref

                region_annotation = db.query_sequence(ggroup_id)
                # logger.info("REGION_ANNOTATION: %s (%s)", str(region_annotation), ggroup_id)
                if region_annotation is not None:
                    _, _, region_annotation = region_annotation
                    category_features = dict(region_annotation).get(str(category.id))
                    if category_features is not None:
                        category_counts[0] += counts  # category row
                        for cf in category_features:
                            category_counts[category_index.get(int(cf))] += counts

                elif it == 0:
                    self.unannotated_counts += counts[:4]
            
            count_sums = counter.counts[1:].sum(axis=0)

            uniq_scaling_factor = (count_sums[0] / count_sums[2], 1.0)[count_sums[2] == 0]
            ambig_scaling_factor = (count_sums[1] / count_sums[3], 1.0)[count_sums[3] == 0]

            logger.info(
                "GCA:: %s CATEGORY COUNTS: uraw=%s unorm=%s araw=%s anorm=%s => SF: %s %s",
                category.name,
                count_sums[0], count_sums[2], count_sums[1], count_sums[3],
                uniq_scaling_factor, ambig_scaling_factor,            
            )

            yield category.name, category_counts, category_index, category_names, uniq_scaling_factor, ambig_scaling_factor

                




    def annotate(self, refmgr, db: AnnotationDatabaseManager, counter: AlignmentCounter, gene_group_db=False):
        """ Annotate a set of gene counts with functional annotations. """
        # self.total_gene_counts, u_sf, a_sf = counter.generate_gene_count_matrix(refmgr)
        # logger.info("TOTAL_GENE_COUNTS = %s", self.total_gene_counts)

        # writer.write_gene_counts(
        #     counter,
        #     refmgr,
        #     u_sf, a_sf,
        #     gene_group_db=gene_group_db,
        # )

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
                gene_id, ggroup_id = rid, rid
            else:
                ref, _ = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)
                gene_id, ggroup_id = ref, ref

            region_annotation = db.query_sequence(ggroup_id)
            if region_annotation is not None:
                _, _, region_annotation = region_annotation
                # logger.info(
                #     "GCAnnotator: Distributing counts of Gene %s (group=%s) %s %s",
                #     gene_id, ggroup_id, counts[0], counts[2],
                # )
                self.distribute_feature_counts(counts, region_annotation)

            else:
                # logger.info("GCAnnotator: Gene %s (group=%s) has no information in database.", gene_id, ggroup_id)
                self.unannotated_counts += counts[:4]

        self.calculate_scaling_factors()

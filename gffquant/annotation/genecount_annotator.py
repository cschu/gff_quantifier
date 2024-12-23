""" module docstring """
import logging

from .count_annotator import CountAnnotator
from .count_writer import CountWriter
from ..counters import AlignmentCounter


logger = logging.getLogger(__name__)


class GeneCountAnnotator(CountAnnotator):
    """ CountAnnotator subclass for gene-based counting. """

    def __init__(self, strand_specific, report_scaling_factors=True):
        """ __init__() """
        CountAnnotator.__init__(self, strand_specific, report_scaling_factors=report_scaling_factors)

    def annotate(self, refmgr, db, counter: AlignmentCounter, gene_group_db=False):
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
            
            ref, _ = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)

            if gene_group_db:
                # ref_tokens = ref.split(".")
                # gene_id, ggroup_id = ".".join(ref_tokens[:-1]), ref_tokens[-1]
                gene_id, ggroup_id = rid, rid
            else:
                gene_id, ggroup_id = ref, ref

            region_annotation = db.query_sequence(ggroup_id)
            if region_annotation is not None:
                _, _, region_annotation = region_annotation
                logger.info(
                    "GCAnnotator: Distributing counts of Gene %s (group=%s) %s %s",
                    gene_id, ggroup_id, counts[0], counts[2],
                )
                self.distribute_feature_counts(counts, region_annotation)

            else:
                logger.info("GCAnnotator: Gene %s (group=%s) has no information in database.", gene_id, ggroup_id)
                self.unannotated_counts += counts[:4]

        self.calculate_scaling_factors()

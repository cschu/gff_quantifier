import logging 

import numpy as np

from .count_annotator import CountAnnotator
from ..counters import CountManager, AlignmentCounter


logger = logging.getLogger(__name__)


class GeneCountAnnotator(CountAnnotator):
	""" CountAnnotator subclass for gene-based counting. """

	def __init__(self, strand_specific, report_scaling_factors=True):
		CountAnnotator.__init__(self, strand_specific, report_scaling_factors=report_scaling_factors)

	def annotate(self, refmgr, db, counter: AlignmentCounter, gene_group_db=False):
		self.total_gene_counts = counter.transform(refmgr)  # count_manager.transform_counts(refmgr)
		logger.info("TOTAL_GENE_COUNTS = %s", self.total_gene_counts)
		# self.total_counts = self.total_gene_counts  # ?

		for rid in counter.get_all_regions():
			counts = counter.get_counts(rid)
			ref, _ = refmgr.get(rid[0] if isinstance(rid, tuple) else rid)

			if gene_group_db:
				ref_tokens = ref.split(".")
				gene_id, ggroup_id = ".".join(ref_tokens[:-1]), ref_tokens[-1]
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


			

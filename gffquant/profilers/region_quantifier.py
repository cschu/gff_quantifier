# pylint: disable=C0103

""" module docstring """

import logging

from ..db.annotation_db import AnnotationDatabaseManager
from .feature_quantifier import FeatureQuantifier

from .. import __tool__, DistributionMode, RunMode

logger = logging.getLogger(__name__)


class RegionQuantifier(FeatureQuantifier):
    # pylint: disable=R0913,R0801
    def __init__(
        self,
        db=None,
        out_prefix=__tool__,
        distribution_mode=DistributionMode.ONE_OVER_N,
        strand_specific=False,
        paired_end_count=1,
        run_mode=RunMode.DOMAIN,
        calculate_coverage=False,
    ):
        FeatureQuantifier.__init__(
            self,
            db=db,
            out_prefix=out_prefix,
            distribution_mode=distribution_mode,
            strand_specific=strand_specific,
            run_mode=run_mode,
            paired_end_count=paired_end_count,
            calculate_coverage=calculate_coverage,
        )

        self.adm = AnnotationDatabaseManager.from_db(self.db)

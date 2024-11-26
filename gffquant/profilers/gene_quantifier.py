# pylint: disable=C0103,R0917

""" module docstring """

import logging

from .feature_quantifier import FeatureQuantifier
from .. import __tool__, DistributionMode, RunMode

logger = logging.getLogger(__name__)


class GeneQuantifier(FeatureQuantifier):
    # pylint: disable=R0913,W0511,R0801
    def __init__(
        self,
        db=None,
        out_prefix=__tool__,
        distribution_mode=DistributionMode.ONE_OVER_N,
        strand_specific=False,
        paired_end_count=1,
        debug=False,
    ):
        FeatureQuantifier.__init__(
            self,
            db=db,
            out_prefix=out_prefix,
            distribution_mode=distribution_mode,
            strand_specific=strand_specific,
            run_mode=RunMode.GENE,
            paired_end_count=paired_end_count,
            debug=debug,
        )

# pylint: disable=C0103

""" module docstring """

import logging

from .feature_quantifier import FeatureQuantifier
from .. import __tool__, RunMode

logger = logging.getLogger(__name__)


class GeneQuantifier(FeatureQuantifier):
    # pylint: disable=R0913,W0511,R0801
    def __init__(
        self,
        db=None,
        out_prefix=__tool__,
        ambig_mode="uniq_only",
        strand_specific=False,
        paired_end_count=1,
    ):
        FeatureQuantifier.__init__(
            self,
            db=db,
            out_prefix=out_prefix,
            ambig_mode=ambig_mode,
            strand_specific=strand_specific,
            run_mode=RunMode.GENE,
            paired_end_count=paired_end_count,
        )

# pylint: disable=C0103

""" module docstring """

import logging

from ..db.annotation_db import AnnotationDatabaseManager
from .feature_quantifier import FeatureQuantifier

from .. import __tool__

logger = logging.getLogger(__name__)


class RegionQuantifier(FeatureQuantifier):
    # pylint: disable=R0913,R0801
    def __init__(
        self,
        db=None,
        out_prefix=__tool__,
        ambig_mode="uniq_only",
        strand_specific=False,
        paired_end_count=1,
        reference_type="genome",
    ):
        FeatureQuantifier.__init__(
            self,
            db=db,
            out_prefix=out_prefix,
            ambig_mode=ambig_mode,
            strand_specific=strand_specific,
            reference_type=reference_type,
            paired_end_count=paired_end_count,
        )

        self.adm = AnnotationDatabaseManager.from_db(self.db)

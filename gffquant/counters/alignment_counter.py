import os
from collections import Counter

class AlignmentCounter(Counter):
	COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled"]

	@staticmethod
	def normalise_counts(counts, feature_len, scaling_factor):
		'''Returns raw, length-normalised, and scaled feature counts.'''
		normalised = counts / feature_len
		scaled = normalised * scaling_factor
		return counts, normalised, scaled
	
	def __init__(self, distribution_mode="uniq_only"):
		Counter.__init__(self)
		self.distribution_mode = distribution_mode

	def dump(self, prefix):
		with open(f"{prefix}.{self.__class__.__name__}.tsv", "wt") as _out:
			for k, v in self.items():
				print(k, v, sep="\t", file=_out)


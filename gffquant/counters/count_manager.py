from .region_counter import RegionCounter
from .seq_counter import UniqueSeqCounter, AmbiguousSeqCounter

class CountManager:
	def __init__(
		self,
		distribution_mode="1overN",
		region_counts=True,
		strandedness_required=False
	):
		self.distribution_mode = distribution_mode
		self.strandedness_required = strandedness_required

		self.uniq_seqcounts = UniqueSeqCounter(
			strandedness_required=strandedness_required
		)

		self.ambig_seqcounts = None

		if distribution_mode not in ("uniq_only", "primary_only"):
			self.ambig_seqcounts = AmbiguousSeqCounter(
				strandedness_required=strandedness_required,
				distribution_mode=distribution_mode
			)

		self.uniq_regioncounts, self.ambig_regioncounts = None, None

		if region_counts:
			self.uniq_regioncounts = RegionCounter(
				strandedness_required=strandedness_required
			)
			if distribution_mode in ("1overN",):	
				self.ambig_regioncounts = RegionCounter(
					strandedness_required=strandedness_required,
					distribution_mode=distribution_mode
				)
	
	@staticmethod
	def _windup_stream(stream):
		for obj in stream:
			...

	def update_unique_counts(self, count_stream):
		stream = self.uniq_seqcounts.update_counts(count_stream)
		if self.uniq_regioncounts is not None:
			stream = self.uniq_regioncounts.update_counts(stream)
		
		CountManager._windup_stream(stream)

	
	def update_ambiguous_counts(self, count_stream):
		stream = self.ambig_seqcounts.update_counts(
			count_stream,
			uniq_counts=self.uniq_seqcounts
		)
		if self.ambig_regioncounts is not None:
			stream = self.ambig_regioncounts.update_counts(stream)

		CountManager._windup_stream(stream)

	def dump_raw_counters(self, prefix):
		if self.uniq_seqcounts is not None:
			self.uniq_seqcounts.dump(prefix)
		if self.ambig_seqcounts: is not None:
			self.ambig_seqcounts.dump(prefix)
		if self.uniq_regioncounts is not None:
			self.uniq_regioncounts.dump(prefix)
		if self.ambig_regioncounts is not None:
			self.ambig_regioncounts.dump(prefix)

		
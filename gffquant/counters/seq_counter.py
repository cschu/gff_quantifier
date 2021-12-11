from .alignment_counter import AlignmentCounter

class UniqueSeqCounter(AlignmentCounter):
	def __init__(
		self,
		strandedness_required=False
	):
		AlignmentCounter.__init__(self)
		self.strandedness_required = strandedness_required

	def get_counts(self, seq_ids):
		if self.strandedness_required:
			return sum(
				self[(seq_id, strand)]
				for seq_id in seq_ids
				for strand in (True, False)
			)
		return sum(self[seq_id] for seq_id in seq_ids)

	def update_counts(self, count_stream):
		increment = 1
		for counts, aln_count, unaligned in count_stream:
			for rid, hits in counts.items():
				updates = []
				if self.strandedness_required:
					strands = {
						strand
						for _, _, strand, _, _ in hits					
					}
					updates += (((rid, strand), increment) for strand in strands)
					
				else:
					updates.append((rid, increment))

				self.update(updates)

			yield counts, aln_count, unaligned


class AmbiguousSeqCounter(AlignmentCounter):
	def __init__(
		self,
		strandedness_required=False,
		distribution_mode="1overN"
	):
		AlignmentCounter.__init__(self, distribution_mode=distribution_mode)
		self.strandedness_required = strandedness_required

	def update_counts(self, count_stream, uniq_counts=None):
		if not uniq_counts:
			raise ValueError("No unique counts supplied.")
		
		def get_increment(uniq_counts, n, n_aln, distribution_mode):
			"""
			currently:
			use dist1 for seq count distribution if there are relevant unique counts
			"""
			if distribution_mode == "all1":
				return 1
			return uniq_counts / n * n_aln if uniq_counts and n else 1 / n_aln
	
		for counts, aln_count, unaligned in count_stream:
			n_total = uniq_counts.get_counts(counts.keys())

			for rid, hits in counts.items():
				updates = []
				if self.strandedness_required:
					strands = {
						strand
						for _, _, strand, _, _ in hits					
					}
					updates += (
						(
							(rid, strand),
							get_increment(
								uniq_counts[(rid, strand)], 
								n_total, aln_count, self.distribution_mode)
						)
						for strand in strands
					)
				else:
					updates.append(
						(rid, get_increment(
							uniq_counts[rid], n_total, aln_count,
							self.distribution_mode))
					)

				self.update(updates)

			yield counts, aln_count, unaligned
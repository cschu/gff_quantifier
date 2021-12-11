from ..bamreader import BamFile, SamFlags

class PairedEndAlignmentCache(dict):
	def __init__(self, ambig_alignments=False):
		dict.__init__(self)
		self.ambig_alignments = ambig_alignments

	def empty_cache(self):
		for qname, alignment in self.items():
			n_aln = len(alignment)

			if n_aln == 1:
				# unpaired
				start, end, flag = alignment[0][1:]
				rev_strand = SamFlags.is_reverse_strand(flag)
			elif n_aln == 2:
				# paired
				start, end = BamFile.calculate_fragment_borders(*alignment[0][1:-1], *alignment[1][1:-1])
				rev_strand = None  # TODO: add strand-specific handling by RNAseq protocol
			else:
				raise ValueError(f"{n_aln} primary alignments detected for read {qname}.")

			yield alignment[0][0], start, end, rev_strand

		self.clear()

	def process_alignment(self, alignment):
		# need to keep track of read pairs to avoid dual counts

		start, end, rev_strand = None, None, None

		# check if the mate has already been seen
		mates = self.setdefault(alignment.qname, list())
		if mates:
			if alignment.rnext != mates[0][0]:
				raise ValueError("Alignment ${alignment.qname} seems to be corrupted: {str(alignment)} {str(mates[0])}.")

			# if mate has been seen, calculate the total fragment size and remove the pair from the cache
			start, end = BamFile.calculate_fragment_borders(alignment.start, alignment.end, *mates[0][1:-1])
			del self[alignment.qname]
		else:
			# otherwise cache the first encountered mate and advance to the next read
			mates.append((alignment.rid, alignment.start, alignment.end, alignment.flag))

		return start, end, rev_strand

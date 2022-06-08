import pysam

from gffquant.bamreader import BamAlignment


class AlignmentProcessor:
	def __init__(self, aln_source="-", aln_type="bam"):
		aln_type = aln_type.lower()
		assert aln_type in ("bam", "sam", "cram")

		self.used_refs = {}
		self.aln_stream = pysam.AlignmentFile(aln_source, "rb" if aln_type == "bam" else "r")

	def get_reference(self, rid):
		return self.used_refs.get(rid, (None, None))

	def get_alignments(self, min_identity=0.97, min_seqlen=45, allow_multiple=True, allow_unique=True, filter_flags=0, required_flags=0):
		with self.aln_stream:
			for aln in self.aln_stream:

				aln = BamAlignment(
					aln.qname,
					aln.flag,
					aln.reference_id,
					aln.pos,
					aln.mapq,
					[(y, x) for x, y in aln.cigar],
					aln.rnext,
					aln.pnext,
					aln.tlen,
					aln.alen,
					dict(aln.tags)
				)

				if aln.flag & filter_flags:
					continue

				if aln.flag & required_flags != required_flags:
					continue

				if (aln.is_ambiguous() and not allow_multiple) or (not aln.is_ambiguous() and not allow_unique):
					continue

				if aln.len_seq < min_seqlen:
					continue

				seqid = 1 - aln.tags.get("NM", 0) / aln.len_seq
				if seqid < min_identity:
					continue

				rname = self.aln_stream.get_reference_name(aln.rid)
				self.used_refs[aln.rid] = rname, self.aln_stream.get_reference_length(rname)

				yield aln

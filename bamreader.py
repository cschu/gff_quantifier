import sys
import time
import gzip
import struct

from collections import Counter

class BamAlignment:
	CIGAR_OPS = "MIDNSHP=X"
	@staticmethod
	def parse_cigar(cigar_ops):
		# op_len << 4 | op
		return [(c >> 4, c & 0xf) for c in cigar_ops]
	@staticmethod
	def show_cigar(cigar):
		return "".join(["{oplen}{op}".format(oplen=c[0], op=BamAlignment.CIGAR_OPS[c[1]]) for c in cigar])
	@staticmethod
	def calculate_coordinates(start, cigar):
		# MIDNSHP=X
		# 012345678; 02378 consume reference: 0000 0010 0011 0111 1000
		REF_CONSUMERS = {0, 2, 3, 7, 8}
		return start + sum(oplen for oplen, op in cigar if op in REF_CONSUMERS)
	def is_primary(self):
		return not self.flag & 0x100
	def is_supplementary(self):
		return self.flag & 0x800
	def is_paired(self):
		return self.flag & 0x1
	def __init__(self, qname=None, flag=None, rid=None, pos=None,
				 mapq=None, cigar=None, rnext=None, pnext=None,
				 tlen=None, len_seq=None, tags=None):
		self.qname = qname
		self.flag = flag
		self.rid = rid
		self.cigar = BamAlignment.parse_cigar(cigar)
		self.start = pos
		self.end = BamAlignment.calculate_coordinates(self.start, self.cigar)
		self.mapq = mapq
		self.rnext = rnext
		self.pnext = pnext
		self.tlen = tlen
		self.len_seq = len_seq
		self.tags = tags
	def __str__(self):
		return "{rid}:{rstart}-{rend} ({cigar};{flag};{mapq};{tlen}) {rnext}:{pnext}".format(
			rid=self.rid, rstart=self.start, rend=self.end, cigar=BamAlignment.show_cigar(self.cigar), flag=self.flag,
			mapq=self.mapq, tlen=self.tlen, rnext=self.rnext, pnext=self.pnext
		 )

class BamFile:
	def __init__(self, fn):
		self._references = list()
		self._file = gzip.open(fn, "rb")
		self._fpos = 0
		self._data_block_start = 0
		self._read_header()

	def _read_header(self):
		try:
			magic = "".join(map(bytes.decode, struct.unpack("cccc", self._file.read(4))))
		except:
			raise ValueError("Could not infer file type.")
		if not magic.startswith("BAM"):
			raise ValueError("Not a valid bam file.")
		len_header = struct.unpack("I", self._file.read(4))[0]
		if len_header:
			header_str = struct.unpack("c" * len_header, self._file.read(len_header))
			# print(header_str)
		n_ref = struct.unpack("I", self._file.read(4))[0]
		self._fpos += 4 + 4 + len_header + 4
		for i in range(n_ref):
			len_rname = struct.unpack("I", self._file.read(4))[0]
			rname = "".join(map(bytes.decode, struct.unpack("c" * len_rname, self._file.read(len_rname))[:-1]))
			len_ref = struct.unpack("I", self._file.read(4))[0]
			self._references.append((rname, len_ref))
			self._fpos += 4 + len_rname + 4
		self._data_block_start = self._fpos

	def get_reference(self, rid):
		return self._references[rid]
	def get_position(self):
		return self._fpos
	def rewind(self):
		self._fpos = self._data_block_start
		self._file.seek(self._data_block_start)
	def seek(self, pos):
		if pos < self._data_block_start:
			raise ValueError("Position {pos} seeks back into the header (data_start={data_start}).".format(pos=pos, data_start=self._data_block_start))
		self._file.seek(pos)
		self._fpos = pos
	def get_alignments(self, required_flags=None, disallowed_flags=None):
		while True:
			try:
				aln_size = struct.unpack("I", self._file.read(4))[0]
			except struct.error:
				break

			rid, pos, len_rname, mapq, bai_bin, n_cigar_ops, flag, len_seq, next_rid, next_pos, tlen = struct.unpack(
				"iiBBHHHIiii",
				self._file.read(32)
			)

			qname = self._file.read(len_rname) # qname
			# qname = struct.unpack("{size}s".format(size=len_rname), self._file.read(len_rname))[0].decode().rstrip("\x00") # qname
			cigar = struct.unpack("I" * n_cigar_ops, self._file.read(4 * n_cigar_ops))
			self._file.read((len_seq + 1) // 2) # encoded read sequence
			self._file.read(len_seq) # quals

			total_read = 32 + len_rname + 4 * n_cigar_ops + (len_seq + 1) // 2 + len_seq
			tags = self._file.read(aln_size - total_read) # get the tags
			self._fpos += 4 + aln_size
			# print(*map(bytes.decode, struct.unpack("c"*len(tags), tags)))

			# yield rid, self._references[rid][0], pos, len_seq
			flag_check = (required_flags is None or (flag & required_flags)) and (disallowed_flags is None or not (flag & disallowed_flags)) 
			if flag_check:
				yield BamAlignment(qname, flag, rid, pos, mapq, cigar, next_rid, next_pos, tlen, len_seq, tags)
			#yield rid, self._references[rid][0], pos, len_seq, qname, cigar, flag, next_rid, self._references[next_rid][0], next_pos, tlen, tags


	def get_multimappers(self):
		'''Returns information (readname, number of ambiguous alignments, genomic coordinates) on all multimapping reads in the bamfile.'''

		multimappers = dict()

		t0 = time.time()
		# First pass: get the ids of all reads that have secondary alignments
		if self._fpos != self._data_block_start:
			self.rewind()
		reads_with_secaln = Counter(
			aln.qname for aln in self.get_alignments(
				required_flags=0x100, disallowed_flags=0x800
			)
		)
		# Second pass: collect information for each multimapping read
		self.rewind()
		for aln_count, aln in enumerate(self.get_alignments(disallowed_flags=0x900), start=1):
			if aln.qname in reads_with_secaln:
				multimappers.setdefault(aln.qname, [reads_with_secaln[aln.qname], set()])[1].add(
					(aln.rid, aln.start, aln.end)
				)
		self.rewind()
		t1 = time.time()
		print("Collected information for {n_multimappers} secondary alignments ({size} bytes) in {n_seconds:.3f}s. (total: {n_alignments})".format(
			size=sys.getsizeof(multimappers), n_multimappers=len(multimappers), n_seconds=t1-t0, n_alignments=len(multimappers) + aln_count), flush=True,
		)
		missing = len(set(reads_with_secaln).difference(multimappers))
		if missing:
			print("{n_missing} secondary alignments don't have primary alignment in file.".format(n_missing=missing), flush=True)

		return multimappers

	@staticmethod
	def calculate_fragment_borders(start1, end1, start2, end2):
		'''
		Returns the minimum and maximum coordinates of a read pair.
		'''
		coords = sorted((start1, end1, start2, end2))
		return coords[0], coords[-1]


if __name__ == "__main__":
	bam = BamFile(sys.argv[1])
	for aln in bam.get_alignments():
		print(aln)
	#print(next(bam.get_alignments()))

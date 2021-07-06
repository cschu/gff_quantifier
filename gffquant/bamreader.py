import sys
import time
import gzip
import struct
import hashlib
import re

# https://stackoverflow.com/questions/22216076/unicodedecodeerror-utf8-codec-cant-decode-byte-0xa5-in-position-0-invalid-s

class SamFlags:
	PAIRED = 0x1
	MATE_UNMAPPED = 0x8
	REVERSE = 0x10
	FIRST_IN_PAIR = 0x40
	SECOND_IN_PAIR = 0x80
	SECONDARY_ALIGNMENT = 0x100
	SUPPLEMENTARY_ALIGNMENT = 0x800

	@staticmethod
	def is_reverse_strand(flag):
		return bool(flag & SamFlags.REVERSE)

class CigarOps:
	CIGAR_OPS = "MIDNSHP=X"
	# MIDNSHP=X
	# 012345678; 02378 consume reference: 0000 0010 0011 0111 1000
	REF_CONSUMERS = {0, 2, 3, 7, 8}

	@staticmethod
	def parse_cigar(cigar_ops):
		# op_len << 4 | op
		return [(c >> 4, c & 0xf) for c in cigar_ops]
	@staticmethod
	def show_cigar(cigar):
		return "".join(["{oplen}{op}".format(oplen=c[0], op=CigarOps.CIGAR_OPS[c[1]]) for c in cigar])
	@staticmethod
	def calculate_coordinates(start, cigar):
		return start + sum(oplen for oplen, op in cigar if op in CigarOps.REF_CONSUMERS)

class BamAlignment:
	TAG_PARAMS = {
		"A": ("c", 1),
		"c": ("b", 1),  # int8
		"C": ("B", 1),  # uint8
		"s": ("h", 2),  # int16
		"S": ("H", 2),  # uint16
		"i": ("i", 4),  # int32
		"I": ("I", 4),  # uint32
		"f": ("f", 4),  # float
	}

	def is_ambiguous(self):
		return self.mapq == 0
	def is_primary(self):
		return not bool(self.flag & SamFlags.SECONDARY_ALIGNMENT)
	def is_supplementary(self):
		return bool(self.flag & SamFlags.SUPPLEMENTARY_ALIGNMENT)
	def is_unique(self):
		return self.mapq != 0
	def is_reverse(self):
		return bool(self.flag & SamFlags.REVERSE)
	def is_paired(self):
		return bool(self.flag & SamFlags.PAIRED)
	def __init__(self, qname=None, flag=None, rid=None, pos=None,
				 mapq=None, cigar=None, rnext=None, pnext=None,
				 tlen=None, len_seq=None, tags=None):
		self.qname = qname
		self.flag = flag
		self.rid = rid
		self.cigar = CigarOps.parse_cigar(cigar)
		self.start = pos
		self.end = CigarOps.calculate_coordinates(self.start, self.cigar)
		self.mapq = mapq
		self.rnext = rnext
		self.pnext = pnext
		self.tlen = tlen
		self.len_seq = len_seq
		self.tags = tags
	def get_hash(self):
		md5 = hashlib.md5()
		md5.update(self.qname)
		return int(md5.hexdigest(), 16)
	def __str__(self):
		return "{rid}:{rstart}-{rend} ({cigar};{flag};{mapq};{tlen}) {rnext}:{pnext}".format(
			rid=self.rid, rstart=self.start, rend=self.end, cigar=CigarOps.show_cigar(self.cigar),
			flag=self.flag, mapq=self.mapq, tlen=self.tlen, rnext=self.rnext, pnext=self.pnext
		 )

class BamFile:
	def __init__(self, fn, large_header=False):
		self._references = list()
		self._file = gzip.open(fn, "rb")
		self._fpos = 0
		self._data_block_start = 0
		# data structures to deal with large headers
		self._refmap = dict()
		self._reverse_references = dict()

		self._read_header(keep_refs=set() if large_header else None)
		if large_header:
			t0 = time.time()
			print("Screening bam for used reference sequences... ", flush=True, end="")
			present_refs = set(aln.rid for _, aln in self.get_alignments())
			t1 = time.time()
			print(" done. ({}s)".format(t1-t0), flush=True)
			self._file.seek(0)
			self._refmap = {rid: i for i, rid in enumerate(sorted(present_refs))}
			self._read_header(keep_refs=present_refs)
			self._reverse_references = {self.get_reference(rid)[0]: rid for rid in self._refmap}

	def revlookup_reference(self, ref):
		return self._reverse_references.get(ref)

	def _read_references(self, keep_refs=None):
		n_ref = struct.unpack("I", self._file.read(4))[0]
		for i in range(n_ref):
			len_rname = struct.unpack("I", self._file.read(4))[0]
			rname = self._file.read(len_rname)[:-1]
			len_ref = struct.unpack("I", self._file.read(4))[0]
			self._fpos += 8 + len_rname
			if keep_refs is None or i in keep_refs:
				yield rname, len_ref
		return list()

	def _read_header(self, keep_refs=None):
		try:
			magic = "".join(map(bytes.decode, struct.unpack("cccc", self._file.read(4))))
		except Exception as exc:
			raise ValueError("Could not infer file type.") from exc
		if not magic.startswith("BAM"):
			raise ValueError("Not a valid bam file.")
		len_header = struct.unpack("I", self._file.read(4))[0]
		if len_header:
			# since we're quite close to the start, we don't lose much by seeking in bgzip'd file
			# obviously, if the header is to be retained, need to handle differently
			self._file.seek(len_header, 1)
		self._fpos = 12 + len_header
		self._references = list(self._read_references(keep_refs=keep_refs))
		self._data_block_start = self._fpos

	def get_reference(self, rid):
		rid_index = self._refmap.get(rid, rid)
		rname, rlen = self._references[rid_index]
		return rname.decode(), rlen

	def n_references(self):
		return len(self._references)

	def get_refdata(self):
		return {rid: self.get_reference(rid) for rid in self._refmap}

	def get_position(self):
		return self._fpos

	def rewind(self):
		self._fpos = self._data_block_start
		self._file.seek(self._data_block_start)

	def seek(self, pos):
		if pos < self._data_block_start:
			raise ValueError(
				f"Position {pos} seeks back into the header (data_start={self._data_block_start})."
			)
		self._file.seek(pos)
		self._fpos = pos

	def _parse_tags(self, size):
		tags = dict()
		bytes_read = 0

		while bytes_read < size:
			# print(bytes_read, size)
			tag, tag_type = struct.unpack("2sc", self._file.read(3))
			# print(tag, tag_type)
			tag, tag_type = map(bytes.decode, (tag, tag_type))
			tag = "".join(tag)
			params = BamAlignment.TAG_PARAMS.get(tag_type)
			if params is not None:
				fmt, tsize = params
				tags[tag] = struct.unpack(fmt, self._file.read(tsize))
				if tsize == 1:
					tags[tag] = tags[tag][0]
			elif tag_type in ("Z", "H"):
				tsize = 0 # we're directly updating bytes_read
				tag_str = list()
				while bytes_read < size:
					_byte = self._file.read(1)
					try:
						_byte = struct.unpack("c", _byte)
					except TypeError:
						break
					_byte = _byte[0].decode()
					if _byte == "\x00":
						bytes_read += 1
						break
					tag_str.append(_byte)
					bytes_read += 1
				tags[tag] = "".join(tag_str)
			elif tag_type == "B":
				tag_type, array_size = struct.unpack("cI", self._file.read(5))
				params = BamAlignment.TAG_PARAMS.get(tag_type.decode())
				tsize = array_size * params[1]
				tags[tag] = struct.unpack(f"{array_size}{params[0]}", self._file.read(tsize))
				tsize += 5
			else:
				raise ValueError(f"Cannot work with tag type {tag}:{tag_type}")
			bytes_read += tsize + 3

		return tags


	def get_alignments(
		self, required_flags=None, disallowed_flags=None, allow_unique=True, allow_multiple=True,
		min_identity=None, min_seqlen=None
	):
		aln_count = 1
		if not allow_multiple and not allow_unique:
			raise ValueError("Either allow_multiple or allow_unique needs to be True.")

		while True:
			try:
				aln_size = struct.unpack("I", self._file.read(4))[0]
			except struct.error:
				break

			rid, pos, len_rname, mapq, _, n_cigarops, flag, len_seq, next_rid, next_pos, tlen = struct.unpack(
				"iiBBHHHIiii",
				self._file.read(32)
			)

			# print(*zip(
			#	"rid pos len_rname mapq _ n_cigarops flag len_seq next_rid next_pos tlen".split(" "), 
			#	(rid, pos, len_rname, mapq, _, n_cigarops, flag, len_seq, next_rid, next_pos, tlen)
			# ))

			qname = self._file.read(len_rname)
			# print('qname', qname)
			cigar = struct.unpack("I" * n_cigarops, self._file.read(4 * n_cigarops))
			# print('cigar', cigar)
			self._file.read((len_seq + 1) // 2) # encoded read sequence
			# print('seq', _)
			self._file.read(len_seq) # quals
			# print('quals', _)

			total_read = 32 + len_rname + 4 * n_cigarops + (len_seq + 1) // 2 + len_seq
			self._fpos += 4 + total_read
			tags_size = aln_size - total_read
			tags = self._parse_tags(tags_size)
			self._fpos += tags_size #4 + aln_size

			md_tag = tags.get("MD", "")
			if md_tag:
				len_seq = sum((
					sum(map(lambda x:int(x.group()), re.finditer("[0-9]+", md_tag))),
					sum(map(lambda x:len(x.group()), re.finditer("[A-Z]+", md_tag)))
				))
			else:
				len_seq = None

			if not len_seq:
				print(f"Read {qname} does not have length information. Skipping", file=sys.stderr, flush=True)

			flag_check = all([
				(required_flags is None or (flag & required_flags)),
				(disallowed_flags is None or not (flag & disallowed_flags))
			])

			is_unique = mapq != 0
			atype_check = any([
				(is_unique and allow_unique),
				(not is_unique and allow_multiple)
			])

			filter_check = all([
				len_seq,
				(min_seqlen is None or len_seq >= min_seqlen),
				(min_identity is None or 1 - (tags.get("NM", 0) / len_seq) >= min_identity)
			])

			if all([flag_check, atype_check, filter_check]):
				yield (
					aln_count,
					BamAlignment(qname, flag, rid, pos, mapq, cigar, next_rid, next_pos, tlen, len_seq, tags)
				)
			aln_count += 1

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

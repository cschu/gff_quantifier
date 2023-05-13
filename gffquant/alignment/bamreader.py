# pylint: disable=C0103,C0301
# pylint: disable=R0902,R0903,R0913,R0914

""" module docstring """

import contextlib
import sys
import time
import gzip
import struct
import hashlib
import re


# https://stackoverflow.com/questions/22216076/unicodedecodeerror-utf8-codec-cant-decode-byte-0xa5-in-position-0-invalid-s


class SamFlags:
    PAIRED = 0x1
    PROPERLY_PAIRED = 0x2
    UNMAPPED = 0x4
    MATE_UNMAPPED = 0x8
    REVERSE = 0x10
    MATE_REVERSE = 0x20
    FIRST_IN_PAIR = 0x40
    SECOND_IN_PAIR = 0x80
    SECONDARY_ALIGNMENT = 0x100
    QUAL_CHECK_FAILURE = 0x200
    PCR_OPTICAL_DUPLICATE = 0x400
    SUPPLEMENTARY_ALIGNMENT = 0x800

    @staticmethod
    def is_reverse_strand(flag):
        return bool(flag & SamFlags.REVERSE)

    @staticmethod
    def is_unmapped(flag):
        return bool(flag & SamFlags.UNMAPPED)


class CigarOps:
    CIGAR_OPS = "MIDNSHP=X"
    #  MIDNSHP=X
    #  012345678; 02378 consume reference: 0000 0010 0011 0111 1000
    REF_CONSUMERS = {0, 2, 3, 7, 8}

    @staticmethod
    def parse_cigar(cigar_ops):
        #  op_len << 4 | op
        return [(c >> 4, c & 0xF) for c in cigar_ops]

    @staticmethod
    def show_cigar(cigar):
        return "".join([f"{c[0]}{CigarOps.CIGAR_OPS[c[1]]}" for c in cigar])

    @staticmethod
    def calculate_coordinates(start, cigar):
        return start + CigarOps.calculate_seqlength(cigar)

    @staticmethod
    def calculate_seqlength(cigar):
        return sum(oplen for oplen, op in cigar if op in CigarOps.REF_CONSUMERS)


class BamAlignment:
    TAG_PARAMS = {
        "A": ("c", 1),
        "c": ("b", 1),  # int8
        "C": ("B", 1),  # uint8
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

    def is_first(self):
        return self.flag & SamFlags.FIRST_IN_PAIR == SamFlags.FIRST_IN_PAIR

    def is_second(self):
        return self.flag & SamFlags.SECOND_IN_PAIR == SamFlags.SECOND_IN_PAIR

    @classmethod
    def from_pysam_alignment(cls, pysam_aln):
        return cls(
            pysam_aln.qname,
            pysam_aln.flag,
            pysam_aln.reference_id,
            pysam_aln.pos,
            pysam_aln.mapq,
            [(y, x) for x, y in pysam_aln.cigar],
            pysam_aln.rnext,
            pysam_aln.pnext,
            pysam_aln.tlen,
            pysam_aln.alen,
            dict(pysam_aln.tags)
        )

    def __init__(
        self,
        qname=None,
        flag=None,
        rid=None,
        pos=None,
        mapq=None,
        cigar=None,
        rnext=None,
        pnext=None,
        tlen=None,
        len_seq=None,
        tags=None,
    ):
        self.qname = qname
        self.flag = flag
        self.rid = rid
        self.cigar = cigar
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
        return (
            f"{self.rid}:{self.start}-{self.end}"
            f"({CigarOps.show_cigar(self.cigar)};{self.flag};{self.mapq};{self.tlen})"
            f"{self.rnext}:{self.pnext}"
        )

    def shorten(self):
        return (self.rid, self.start, self.end, self.is_reverse())


class BamBuffer:
    def __init__(self, fn, size=1000000):
        self._file = gzip.open(fn, "rb")
        self._buffer_size = size
        self._buffer = None
        self._fpos = 0
        self._bp = 0
        self._fill()

        self._data_block_start = None

    def _fill(self):
        self._buffer = self._file.read(self._buffer_size)
        self._bp = 0

    def is_exhausted(self):
        return not bool(self._file) and self.is_empty()

    def is_empty(self):
        return self.bytes_in_buffer() == 0

    def bytes_in_buffer(self):
        return len(self._buffer) - self._bp

    def set_data_block_start(self):
        if self._data_block_start is None:
            self._data_block_start = self._fpos
        elif self._data_block_start != self._fpos:
            raise ValueError(
                f"Data block start was already set to {self._data_block_start}."
                f"self._fpos={self._fpos}"
            )

    def skip(self, size):
        _ = self.read(size, skip=True)

    def rewind(self, full=False):
        if full:
            self._file.seek(0)
            self._fpos = 0
        else:
            self._fpos = self._data_block_start
            self._file.seek(self._data_block_start)
        self._fill()

    def read(self, nbytes, skip=False):
        extra_bytes_required = nbytes - self.bytes_in_buffer()
        _bytes = None
        if extra_bytes_required > 0:
            if skip:
                self._file.seek(extra_bytes_required, 1)
                self._fill()
            else:
                _bytes = self._buffer[self._bp:]  # take the remaining bytes from buffer
                _bytes += self._file.read(extra_bytes_required)  # obtain required bytes
                self._fill()  # refill the buffer
        else:
            _bytes = self._buffer[self._bp:self._bp + nbytes]
            self._bp += nbytes

        self._fpos += nbytes
        return _bytes


class BamFile:
    def __init__(self, fn, large_header=False, buffer_size=1000000, ambig_bookkeeper=None):
        self._references = []
        self._file = BamBuffer(fn, size=buffer_size)
        # data structures to deal with large headers
        self._refmap = {}
        self._reverse_references = {}

        self._read_header(keep_refs=set() if large_header else None)
        do_ambig_bookkeeping = ambig_bookkeeper is not None and not isinstance(ambig_bookkeeper, contextlib.nullcontext)

        if large_header or do_ambig_bookkeeping:
            t0 = time.time()
            print("Screening bam for used reference sequences... ", flush=True, end="")
            present_refs = set()
            for _, aln in self.get_alignments(parse_tags=False, reference_screening=True):
                present_refs.add(aln.rid)
                if do_ambig_bookkeeping and aln.is_ambiguous():
                    ambig_bookkeeper.register_alignment(aln)

            # present_refs = set(
            #     aln.rid
            #     for _, aln in self.get_alignments(
            #         parse_tags=False, reference_screening=True
            #     )
            # )
            t1 = time.time()
            print(f" done. ({t1-t0}s)", flush=True)

            self._file.rewind(full=True)
            self._refmap = {rid: i for i, rid in enumerate(sorted(present_refs))}
            self._read_header(keep_refs=present_refs)
            self._reverse_references = {
                self.get_reference(rid)[0]: rid for rid in self._refmap
            }

    def revlookup_reference(self, ref):
        return self._reverse_references.get(ref)

    def _read_references(self, keep_refs=None):
        n_ref = struct.unpack("I", self._file.read(4))[0]

        for i in range(n_ref):
            len_rname = struct.unpack("I", self._file.read(4))[0]
            rname = self._file.read(len_rname)[:-1]
            len_ref = struct.unpack("I", self._file.read(4))[0]
            if keep_refs is None or i in keep_refs:
                yield rname, len_ref

        return []

    def _read_header(self, keep_refs=None):
        try:
            magic = "".join(
                map(bytes.decode, struct.unpack("cccc", self._file.read(4)))
            )
        except Exception as exc:
            raise ValueError("Could not infer file type.") from exc
        if not magic.startswith("BAM"):
            raise ValueError("Not a valid bam file.")
        len_header = struct.unpack("I", self._file.read(4))[0]
        if len_header:
            # since we're quite close to the start, we don't lose much by seeking in bgzip'd file
            # obviously, if the header is to be retained, need to handle differently
            # self._file.seek(len_header, 1)
            self._file.skip(len_header)

        self._references = list(self._read_references(keep_refs=keep_refs))
        self._file.set_data_block_start()

    def get_reference(self, rid):
        rid_index = self._refmap.get(rid, rid)
        rname, rlen = self._references[rid_index]
        return rname.decode(), rlen

    def n_references(self):
        return len(self._references)

    def get_refdata(self):
        return {rid: self.get_reference(rid) for rid in self._refmap}

    @staticmethod
    def _parse_tags(tagdata):
        tags = {}
        p, n = 0, len(tagdata)
        while p < n:
            tag, tag_type = struct.unpack("2sc", tagdata[p:p + 3])
            p += 3
            tag, tag_type = map(bytes.decode, (tag, tag_type))
            tag = "".join(tag)
            params = BamAlignment.TAG_PARAMS.get(tag_type)
            if params is not None:
                fmt, tsize = params
                tags[tag] = struct.unpack(fmt, tagdata[p:p + tsize])
                p += tsize
                if tsize == 1:
                    tags[tag] = tags[tag][0]
            elif tag_type == "B":
                tag_type, array_size = struct.unpack("cI", tagdata[p:p + 5])
                p += 5
                params = BamAlignment.TAG_PARAMS.get(tag_type.decode())
                tsize = array_size * params[1]
                tags[tag] = struct.unpack(
                    f"{array_size}{params[0]}", tagdata[p:p + tsize]
                )
                p += tsize
            elif tag_type in "ZH":
                tag_str, i = [], 0
                for i, _ in enumerate(tagdata[p:]):
                    _byte = tagdata[p + i:p + i + 1]
                    try:
                        _byte = struct.unpack("c", _byte)
                    except TypeError as e:
                        print(f"TypeError: {e}", _byte)
                        break
                    _byte = _byte[0].decode()
                    if _byte == "\x00":
                        p += 1
                        break
                    tag_str.append(_byte)
                tags[tag] = "".join(tag_str)
                p += i
            else:
                raise ValueError(f"Cannot work with tag type {tag}:{tag_type}")

        return tags

    def get_alignments(
        self,
        required_flags=None,
        disallowed_flags=None,
        allow_unique=True,
        allow_multiple=True,
        min_identity=None,
        min_seqlen=None,
        parse_tags=True,
        reference_screening=False,
    ):
        aln_count = 1
        if not allow_multiple and not allow_unique:
            raise ValueError("Either allow_multiple or allow_unique needs to be True.")

        while not self._file.is_exhausted():
            try:
                aln_size = struct.unpack("I", self._file.read(4))[0]
            except struct.error:
                break

            # pylint: disable-next=C0301
            unpacked = struct.unpack("iiBBHHHIiii", self._file.read(32))
            (
                rid,
                pos,
                len_rname,
                mapq,
                _,
                n_cigarops,
                flag,
                len_seq,
                next_rid,
                next_pos,
                tlen,
            ) = unpacked

            qname = self._file.read(len_rname)
            cigar = struct.unpack("I" * n_cigarops, self._file.read(4 * n_cigarops))
            cigar = CigarOps.parse_cigar(cigar)

            if len_seq:
                self._file.read((len_seq + 1) // 2)  # encoded read sequence
                self._file.read(len_seq)  # quals

            total_read = 32 + len_rname + 4 * n_cigarops + (len_seq + 1) // 2 + len_seq
            tags_size = aln_size - total_read
            tag_data = self._file.read(tags_size)
            tags = {}
            if parse_tags or not reference_screening:
                tags = BamFile._parse_tags(tag_data)

            len_seq = CigarOps.calculate_seqlength(cigar)

            if not len_seq:
                md_tag = tags.get("MD")
                if md_tag:
                    len_seq = sum(
                        (
                            sum(
                                map(lambda x: int(x.group()), re.finditer("[0-9]+", md_tag))
                            ),
                            sum(
                                map(lambda x: len(x.group()), re.finditer("[A-Z]+", md_tag))
                            ),
                        )
                    )
                else:
                    len_seq = None

            if parse_tags and not len_seq:
                print(
                    f"Read {qname} does not have length information. Skipping",
                    file=sys.stderr,
                    flush=True,
                )

            # pylint: disable=C0325
            flag_check = all(
                [
                    (required_flags is None or (flag & required_flags)),
                    (disallowed_flags is None or not (flag & disallowed_flags)),
                ]
            )

            is_unique = mapq != 0
            atype_check = any(
                [(is_unique and allow_unique), (not is_unique and allow_multiple)]
            )

            filter_check = all(
                [
                    len_seq,
                    (min_seqlen is None or len_seq >= min_seqlen),
                    (
                        min_identity is None or 1 - (tags.get("NM", 0) / len_seq) >= min_identity
                    ),
                ]
            )

            if reference_screening or all((flag_check, atype_check, filter_check)):
                yield (
                    aln_count,
                    BamAlignment(
                        qname,
                        flag,
                        rid,
                        pos,
                        mapq,
                        cigar,
                        next_rid,
                        next_pos,
                        tlen,
                        len_seq,
                        tags if tags else None,
                    ),
                )
            aln_count += 1

    @staticmethod
    def calculate_fragment_borders(start1, end1, start2, end2):
        """
        Returns the minimum and maximum coordinates of a read pair.
        """
        coords = sorted((start1, end1, start2, end2))
        return coords[0], coords[-1]


if __name__ == "__main__":
    # bam = BamFile(sys.argv[1])
    # for aln in bam.get_alignments():
    #    print(aln)

    bam = BamFile("/Users/cschu/gqdebug/testcase.bam", buffer_size=100, large_header=True)

    L = list(bam.get_alignments())

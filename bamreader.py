import sys
import gzip
import struct

class BamFile:
    def __init__(self, fn):
        self._references = list()
        self._file = gzip.open(fn, "rb")
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
       for i in range(n_ref):
           len_rname = struct.unpack("I", self._file.read(4))[0]
           rname = "".join(map(bytes.decode, struct.unpack("c" * len_rname, self._file.read(len_rname))[:-1]))
           len_ref = struct.unpack("I", self._file.read(4))[0]
           self._references.append((rname, len_ref))
    def get_alignments(self):
        while True:
            try:
                aln_size = struct.unpack("I", self._file.read(4))[0]
            except struct.error:
                break

            rid, pos, len_rname, mapq, bai_bin, n_cigar_ops, flag, len_seq, next_rid, next_pos, tlen = struct.unpack(
                "iiBBHHHIiii",
                self._file.read(32)
            )

            self._file.read(len_rname) # qname
            self._file.read(4 * n_cigar_ops)
            self._file.read((len_seq + 1) // 2) # encoded read sequence
            self._file.read(len_seq) # quals

            total_read = 32 + len_rname + 4 * n_cigar_ops + (len_seq + 1) // 2 + len_seq
            tags = self._file.read(aln_size - total_read) # get the tags
            # print(*map(bytes.decode, struct.unpack("c"*len(tags), tags)))

            yield rid, self._references[rid][0], pos, len_seq



if __name__ == "__main__":
    bam = BamFile(sys.argv[1])
    for aln in bam.get_alignments():
        print(aln)
    #print(next(bam.get_alignments()))

""" CIGAR Ops """


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

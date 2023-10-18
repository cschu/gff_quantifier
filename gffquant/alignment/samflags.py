""" SamFlags (partial) """

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

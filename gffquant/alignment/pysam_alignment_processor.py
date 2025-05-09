# pylint: disable=R0917

""" module docstring """

from contextlib import nullcontext

import pysam

from .bamreader import BamAlignment
from .samflags import SamFlags


class AlignmentProcessor:
    PASS_FILTER = 0
    SEQID_FILTERED = 1
    LENGTH_FILTERED = 2
    TOTAL_READS = 0
    TOTAL_ALIGNED_READS = 1
    TOTAL_PASSED_READS = 2

    def __init__(self, aln_source="-", aln_type="bam"):
        aln_type = aln_type.lower()
        assert aln_type in ("bam", "sam")

        # self.used_refs = {}
        # pylint: disable=E1101
        self.aln_stream = pysam.AlignmentFile(aln_source, "rb" if aln_type == "bam" else "r")
        self.stat_counter = [0, 0, 0]
        self.read_counter = [0, 0, 0]

    # def get_reference(self, rid):
    #     return self.used_refs.get(rid, (None, None))

    def get_alignment_stats(self):
        return self.stat_counter

    def get_alignment_stats_dict(self):
        return dict(
            zip(
                ("pysam_total", "pysam_passed", "pysam_seqid_filt", "pysam_len_filt"),
                [sum(self.stat_counter), ] + self.stat_counter
            )
        )

    @staticmethod
    def get_alignment_stats_str(stat_counter, table=True):
        # pylint: disable=R1705
        if table:
            return "\n".join(
                "\t".join(s)
                for s in zip(
                    ("Total", "Passed", "Seqid", "Length"),
                    (str(v) for v in (sum(stat_counter),) + tuple(stat_counter))
                )
            )
        else:
            return f"Total:{sum(stat_counter)} " + \
                f"Passed filters: {stat_counter[0]} " + \
                f"Filtered(seqid): {stat_counter[1]} " + \
                f"Filtered(length): {stat_counter[2]}"

    def check_alignment(self, pysam_aln, min_seqlen, min_identity):
        if pysam_aln.alen < min_seqlen:
            self.stat_counter[AlignmentProcessor.LENGTH_FILTERED] += 1
            return False

        try:
            mismatches = pysam_aln.get_tag("NM")
        except KeyError:
            mismatches = 0

        seqid = 1 - mismatches / pysam_aln.alen
        if seqid < min_identity:
            self.stat_counter[AlignmentProcessor.SEQID_FILTERED] += 1
            return False

        return True

    # pylint: disable=R0913,R0914,W0613
    def get_alignments(
        self,
        min_identity=0.97,
        min_seqlen=45,
        filter_flags=0,
        required_flags=0,
        filtered_sam=None,
    ):
        last_read, last_aligned_read, last_passed_read = None, None, None

        # pylint: disable=R1732
        filtered_out = open(filtered_sam, "wt", encoding="UTF-8") if filtered_sam else nullcontext()
        # pylint: enable=R1732

        with self.aln_stream, filtered_out:
            for pysam_aln in self.aln_stream:

                read_unmapped = SamFlags.is_unmapped(pysam_aln.flag)

                if last_read is None or pysam_aln.qname != last_read:
                    last_read = pysam_aln.qname
                    self.read_counter[AlignmentProcessor.TOTAL_READS] += 1

                    # if not read_unmapped:
                    #     self.read_counter[AlignmentProcessor.TOTAL_ALIGNED_READS] += 1

                if read_unmapped:
                    continue

                if last_aligned_read is None or pysam_aln.qname != last_aligned_read:
                    last_aligned_read = pysam_aln.qname
                    self.read_counter[AlignmentProcessor.TOTAL_ALIGNED_READS] += 1

                # if read_unmapped:
                #     continue

                if pysam_aln.flag & filter_flags:
                    continue

                if pysam_aln.flag & required_flags != required_flags:
                    continue

                if self.check_alignment(pysam_aln, min_seqlen, min_identity):
                    self.stat_counter[AlignmentProcessor.PASS_FILTER] += 1

                    rname = pysam_aln.reference_name
                    rlength = self.aln_stream.get_reference_length(rname)
                    # self.used_refs[pysam_aln.reference_id] = rname, rlength

                    if last_passed_read is None or pysam_aln.qname != last_passed_read:
                        last_passed_read = pysam_aln.qname
                        self.read_counter[AlignmentProcessor.TOTAL_PASSED_READS] += 1

                    yield BamAlignment.from_pysam_alignment(pysam_aln, rlength)

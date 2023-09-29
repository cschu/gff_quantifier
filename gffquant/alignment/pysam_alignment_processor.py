""" module docstring """

import pysam

from .bamreader import BamAlignment, SamFlags


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

        self.used_refs = {}
        # pylint: disable=E1101
        self.aln_stream = pysam.AlignmentFile(aln_source, "rb" if aln_type == "bam" else "r")
        self.stat_counter = [0, 0, 0]
        self.read_counter = [0, 0, 0]

    def get_reference(self, rid):
        return self.used_refs.get(rid, (None, None))

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

    # pylint: disable=R0913,W0613
    def get_alignments(
        self,
        min_identity=0.97,
        min_seqlen=45,
        allow_multiple=True,
        allow_unique=True,
        filter_flags=0,
        required_flags=0,
        verbose=True,
    ):
        last_read, last_passed_read = None, None

        with self.aln_stream:
            for pysam_aln in self.aln_stream:

                read_unmapped = SamFlags.is_unmapped(pysam_aln.flag)

                if last_read is None or pysam_aln.qname != last_read:
                    last_read = pysam_aln.qname
                    self.read_counter[AlignmentProcessor.TOTAL_READS] += 1

                    if not read_unmapped:
                        self.read_counter[AlignmentProcessor.TOTAL_ALIGNED_READS] += 1

                if read_unmapped:
                    continue

                aln = BamAlignment.from_pysam_alignment(pysam_aln)

                if aln.flag & filter_flags:
                    continue

                if aln.flag & required_flags != required_flags:
                    continue

                if (aln.is_ambiguous() and not allow_multiple):
                    continue

                if (not aln.is_ambiguous() and not allow_unique):
                    continue

                if aln.len_seq < min_seqlen:
                    self.stat_counter[AlignmentProcessor.LENGTH_FILTERED] += 1
                    continue

                seqid = 1 - aln.tags.get("NM", 0) / aln.len_seq
                if seqid < min_identity:
                    self.stat_counter[AlignmentProcessor.SEQID_FILTERED] += 1
                    continue

                rname = self.aln_stream.get_reference_name(aln.rid)
                self.used_refs[aln.rid] = rname, self.aln_stream.get_reference_length(rname)

                self.stat_counter[AlignmentProcessor.PASS_FILTER] += 1

                if last_passed_read is None or pysam_aln.qname != last_passed_read:
                    last_passed_read = pysam_aln.qname
                    self.read_counter[AlignmentProcessor.TOTAL_PASSED_READS] += 1

                yield aln

# pylint: disable=W0223

""" module docstring """

from .alignment_counter import AlignmentCounter


class UniqueSeqCounter(AlignmentCounter):
    def __init__(self, strand_specific=False):
        AlignmentCounter.__init__(self, strand_specific=strand_specific)

    def get_counts(self, seq_ids):
        """
        Given a list of sequence ids, return the total number of reads that mapped to each of those
        sequences

        :param seq_ids: a list of sequence ids to count
        :return: A list of counts for each sequence ID.
        """
        if self.strand_specific:
            return sum(
                self[(seq_id, strand)] for seq_id in seq_ids for strand in (True, False)
            )
        return sum(self[seq_id] for seq_id in seq_ids)

    def update_counts(self, count_stream, increment=1):
        for counts, aln_count, unaligned in count_stream:

            for rid, hits in counts.items():

                if self.strand_specific:
                    strands = tuple(int(strand) for _, _, strand, _, _ in hits)

                    self[(rid, True)] += sum(strands) * increment
                    self[(rid, False)] += (len(hits) - sum(strands)) * increment

                else:
                    self[rid] += len(hits) * increment

            yield counts, aln_count, unaligned


class AmbiguousSeqCounter(AlignmentCounter):
    def __init__(self, strand_specific=False, distribution_mode="1overN"):
        AlignmentCounter.__init__(
            self, distribution_mode=distribution_mode, strand_specific=strand_specific
        )

    def update_counts(self, count_stream):
        def get_increment(n_aln, distribution_mode):
            return 1 / n_aln if distribution_mode == "1overN" else 1

        for counts, aln_count, unaligned in count_stream:

            increment = get_increment(aln_count, self.distribution_mode)

            for rid, hits in counts.items():

                if self.strand_specific:
                    strands = tuple(int(strand) for _, _, strand, _, _ in hits)

                    self[(rid, True)] += sum(strands) * increment
                    self[(rid, False)] += (len(hits) - sum(strands)) * increment

                else:
                    self[rid] += len(hits) * increment
                    print(f"RID:{rid} = {self[rid]} NHITS:{len(hits)}")

            yield counts, aln_count, unaligned

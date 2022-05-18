# pylint: disable=W0105,C0103,W1514,R0902
# pylint: disable=C0301

""" module docstring """

import gzip

from collections import Counter
from dataclasses import dataclass

from gffquant.bamreader import SamFlags

DEBUG = False


@dataclass
class AmbiguousAlignment:
    rid: int
    start: int
    end: int
    cstart: int
    cend: int
    flag: int

    def __hash__(self):
        return hash((self.rid, self.start, self.end, self.cstart, self.cend, self.flag))

    def is_unannotated(self):
        return self.rid == -1

    def is_first(self):
        return self.flag & SamFlags.FIRST_IN_PAIR == SamFlags.FIRST_IN_PAIR

    def is_second(self):
        return self.flag & SamFlags.SECOND_IN_PAIR == SamFlags.SECOND_IN_PAIR

    def is_primary(self):
        return self.flag & SamFlags.SECONDARY_ALIGNMENT == 0


class AmbiguousAlignmentRecordKeeper:
    """
    This class takes care of the specific record keeping for ambiguous alignments.
    This includes:
            - counting of annotated / unannotated reads
            - assignment of integer read ids to save space
            - writing of annotated read information to a temporary file
            (this last one is a painful workaround, the alternative would be to
             split the bam file into unique and ambiguous alignments, with the latter sorted by name)
             The task here is to save memory, which now comes at the cost of a bit of (temporary) disk space.

    The class is used as a contextmanager. It is only instanced when ambiguous alignments require
    special treatment.
    """

    def __init__(self, prefix, db, do_overlap_detection=True):
        self.annotated = set()
        self.unannotated = set()
        self.readids = {}
        self.dumpfile = prefix + ".ambig_tmp.txt.gz"
        # pylint: disable=R1732
        self.ambig_dump = gzip.open(self.dumpfile, "wt")
        self.db = db
        self.do_overlap_detection = do_overlap_detection
        self.aln_counter = Counter()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.ambig_dump.close()

    def register_alignment(self, aln):
        self.aln_counter[(aln.qname, aln.is_second())] += 1

    def get_ambig_alignment_count(self, aln):
        return self.aln_counter.get((aln.qname, aln.is_second()), 1)

    def clear(self):
        self.aln_counter.clear()

    def process_alignment(self, ref, aln, aln_count):
        """
              - obtains/assigns a unique integer (correlating to the read id) to an alignment for
        alignment group identification
              - updates the unannotated/annotated records, depending on whether the alignment could be annotated
              - writes the relevant alignment information to disk if alignment could be annotated
        """
        qname_id = self.readids.setdefault(aln.qname, len(self.readids))
        if self.do_overlap_detection:
            overlaps, coverage = self.db.get_overlaps(ref, aln.start, aln.end)
            if not overlaps:
                self.unannotated.add(qname_id)
            else:
                self.annotated.add(qname_id)
                for ovl, (cstart, cend) in zip(overlaps, coverage):
                    print(
                        qname_id,
                        aln_count,
                        aln.rid,
                        ovl.begin,
                        ovl.end,
                        cstart,
                        cend,
                        aln.flag,
                        file=self.ambig_dump,
                        sep="\t",
                    )
        else:
            print(
                qname_id,
                aln_count,
                aln.rid,
                -1,
                -1,
                -1,
                -1,
                aln.flag,
                file=self.ambig_dump,
                sep="\t",
            )

    def n_unannotated(self):
        """
        returns the number of unannotated reads (all reads that didn't align to any annotated region)
        """
        return len(self.unannotated.difference(self.annotated))


class AmbiguousAlignmentGroup:
    FIRST_MATES, SECOND_MATES = False, True

    def __init__(self, aln):
        self.alignments = [set(), set()]
        self.unannotated = [0, 0]
        self.qname = aln[0]

        self.add_alignment(aln)

    def add_alignment(self, aln):
        aln = AmbiguousAlignment(*aln[2:])

        if aln.is_unannotated():
            self.unannotated[aln.is_first()] += 1
        else:
            self.alignments[aln.is_first()].add(aln)

    def n_align(self):
        return len(self.alignments[0]) + len(self.alignments[1])

    def resolve(self):
        for mate_alignments, ucounts in zip(self.alignments, self.unannotated):
            hits = {}
            for aln in mate_alignments:
                hits.setdefault(aln.rid, set()).add(
                    (aln.start, aln.end, aln.flag, aln.cstart, aln.cend)
                )
            if hits:
                yield hits, len(mate_alignments), ucounts

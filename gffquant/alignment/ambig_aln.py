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
        # aln_data = self.aln_counter.setdefault()
        # qname_id = self.readids.setdefault(aln.qname, len(self.readids))
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
                # print(hits, len(mate_alignments), ucounts)
                yield hits, len(mate_alignments), ucounts

        # hits = {}
        # for aln in chain(*self.alignments):
        #    hits.setdefault(aln.rid, set()).add(
        #        (aln.start, aln.end, aln.flag, aln.cstart, aln.cend)
        #    )

        # return ((hits, self.n_align(), self.unannotated),)


class AmbiguousAlignmentGroup2:
    """
    Represents a group of ambiguous alignments of a single read/read pair.
    bwa -A assigns a primary read (pair) and flags all others as secondary alignments (0x100)
    due to the ngless filtering, we also have the case that the primary was filtered out

    Caveat: paired-end information is currently not considered for ambiguous alignment groups.
    This means that overlapping mates will result in duplicate counts.
    """

    def __init__(self, aln):
        self.primaries = []
        self.alignments = {}
        self.secondaries = {}  # these are the secondary alignments
        self.primary1, self.primary2 = None, None
        self.qname = aln[0]
        self.uniq_alignments = set()  # this is the set of alignments that can be annotated
        self.unannotated = 0
        self.add_alignment(aln)

    def resolve(self):

        hits = {}

        if len(self.primaries) > 2:
            print(
                f"Warning: found more than two primary alignments for read {self.qname}"
            )
            print(*(str(aln) for aln in self.primaries), sep="\n")

        for rid, alignments in self.alignments.items():
            """
            treat_as_individuals = False

            if len(alignments) != 2:
                    treat_as_individuals = True
            else:
                    first = next((aln for aln in alignments if aln.is_first()), None)
                    second = next((aln for aln in alignments if aln.is_second()), None)
                    if first is not None and second is not None:
                            ...
                    else:
                            # in this case, we cannot identify the mate relationship on the same reference sequence
                            # and therefore treat both alignments as individuals
                            treat_as_individuals = True
            """

            for aln in alignments:
                hits.setdefault(rid, set()).add(
                    (aln.start, aln.end, aln.flag, aln.cstart, aln.cend)
                )

        # ugly!!!
        return (
            (
                hits,
                sum(len(alignments) for alignments in self.alignments.values()),
                self.unannotated,
            ),
        )

        # print(self.qname, self.primary1, self.primary2)

        # https://stackoverflow.com/questions/64090762/python-lazy-function-evaluation-in-any-all
        # cannot be evaluated with all
        # if self.primary1 is not None \
        # 	and self.primary2 is not None \
        # 	and self.primary1[:-1] == self.primary2[:-1]:
        # 		self.primary2 = None
        #
        # alignments = set([self.primary1, self.primary2]).union(self.secondaries).difference({None})
        # for rid, start, end, cstart, cend, flag in alignments:
        # 	hits.setdefault(rid, set()).add((start, end, SamFlags.is_reverse_strand(flag), cstart, cend))
        #
        # print("HITS:", hits, self.n_align(), self.unannotated)
        # return ((hits, self.n_align(), self.unannotated),)

    def add_alignment(self, aln):
        print(aln)
        aln = AmbiguousAlignment(*aln[2:])

        if aln.is_unannotated():
            self.unannotated += 1
        else:
            self.alignments.setdefault(aln.rid, []).append(aln)

            if aln.is_primary():
                self.primaries.append(aln)

        """
		processed = False
		if aln.is_unannotated():
			self.unannotated += 1
			processed = True
		elif aln.is_primary():
			if aln.is_first():
				if self.primary1 is None:
					self.primary1 = aln
					processed = True
				else:
					print(f"Warning: found additional first-primary alignment {str(aln)}")

			elif aln.is_second():
				if self.primary2 is None:
					self.primary2 = aln
					processed = True
				else:
					print(f"Warning: found additional second-primary alignment {str(aln)}")

		if not processed:
			self.secondaries.

		rid, start, end, cstart, cend, flag = aln
		"""

        # return None

        # flag = aln[-1]
        # short_aln = tuple(aln[2:])
        # if short_aln[0] == -1:
        # 	self.unannotated += 1
        # elif flag & SamFlags.FIRST_IN_PAIR and self.primary1 is None:
        # 	self.primary1 = short_aln
        # elif flag & SamFlags.SECOND_IN_PAIR and self.primary2 is None:
        # 	self.primary2 = short_aln
        # else:
        # 	self.secondaries.append(short_aln)
        # if short_aln[0] != -1:
        # 	self.uniq_alignments.add(short_aln[:-1])

    def n_align(self):
        return len(
            self.uniq_alignments.union((self.primary1, self.primary2)).difference(
                {None}
            )
        )

    """
	def resolve_old(self):

		hits = {}
		# print(self.qname, self.primary1, self.primary2)

		# https://stackoverflow.com/questions/64090762/python-lazy-function-evaluation-in-any-all
		# cannot be evaluated with all
		if self.primary1 is not None and self.primary2 is not None and self.primary1[:-1] == self.primary2[:-1]:
			self.primary2 = None

		alignments = set([self.primary1, self.primary2]).union(self.secondaries).difference({None})
		for rid, start, end, cstart, cend, flag in alignments:
			hits.setdefault(rid, set()).add((start, end, SamFlags.is_reverse_strand(flag), cstart, cend))

		print("HITS:", hits, self.n_align(), self.unannotated)
		return ((hits, self.n_align(), self.unannotated),)
	"""

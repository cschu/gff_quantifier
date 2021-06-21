from gffquant.bamreader import SamFlags

DEBUG = False


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
		self.readids = dict()
		self.dumpfile = prefix + ".ambig_tmp.txt"
		self.ambig_dump = open(self.dumpfile, "wt")
		self.db = db
		self.do_overlap_detection = do_overlap_detection

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_val, exc_tb):
		self.ambig_dump.close()

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
						qname_id, aln_count, aln.rid, ovl.begin, ovl.end, cstart, cend, aln.flag,
						file=self.ambig_dump, sep="\t"
					)
		else:
			print(
				qname_id, aln_count, aln.rid, -1, -1, -1, -1, aln.flag,
				file=self.ambig_dump, sep="\t"
			)

	def n_unannotated(self):
		"""
		returns the number of unannotated reads (all reads that didn't align to any annotated region)
		"""
		return len(self.unannotated.difference(self.annotated))


class AmbiguousAlignmentGroup:
	"""
	Represents a group of ambiguous alignments of a single read/read pair.
	bwa -A assigns a primary read (pair) and flags all others as secondary alignments (0x100)
	due to the ngless filtering, we also have the case that the primary was filtered out

	Caveat: paired-end information is currently not considered for ambiguous alignment groups.
	This means that overlapping mates will result in duplicate counts.
	"""

	def __init__(self, aln):
		self.secondaries = list()  # these are the secondary alignments
		self.primary1, self.primary2 = None, None
		self.qname = aln[0]
		self.uniq_alignments = set()  # this is the set of alignments that can be annotated
		self.unannotated = 0
		self.add_alignment(aln)

	def add_alignment(self, aln):
		flag = aln[-1]
		short_aln = tuple(aln[2:])
		if short_aln[0] == -1:
			self.unannotated += 1
		elif flag & SamFlags.FIRST_IN_PAIR and self.primary1 is None:
			self.primary1 = short_aln
		elif flag & SamFlags.SECOND_IN_PAIR and self.primary2 is None:
			self.primary2 = short_aln
		else:
			self.secondaries.append(short_aln)
		if short_aln[0] != -1:
			self.uniq_alignments.add(short_aln[:-1])

	def n_align(self):
		return len(self.uniq_alignments.union((self.primary1, self.primary2)).difference({None}))

	def resolve(self, counter, bam, distmode="1overN"):

		hits = dict()
		print(self.qname, self.primary1, self.primary2)

		# https://stackoverflow.com/questions/64090762/python-lazy-function-evaluation-in-any-all 
		# cannot be evaluated with all
		if self.primary1 is not None and self.primary2 is not None and self.primary1[:-1] == self.primary2[:-1]:
			self.primary2 = None

		alignments = set([self.primary1, self.primary2]).union(self.secondaries).difference({None})
		for rid, start, end, cstart, cend, flag in alignments:
			hits.setdefault(rid, set()).add((start, end, cstart, cend, SamFlags.is_reverse_strand(flag)))

		counter.update_ambiguous_counts(
			hits, self.n_align(), self.unannotated, feat_distmode=distmode
		)

import time
import sys
import gzip

from intervaltree import IntervalTree

from bamreader import BamFile

class FeatureQuantifier:
    def _read_gff_index(self, gff_index):
        self.gff_index = dict()
        for line in open(gff_index, "rt"):
            line = line.strip().split("\t")
            self.gff_index.setdefault(line[0], list()).append(list(map(int, line[1:3])))
    def __init__(self, gff_db, gff_index, gff_gzipped=True):
        self.gff_db = gff_db
        gff_open = gzip.open if gff_gzipped else open
        self.gff_data = gff_open(gff_db, "rt")
        self._read_gff_index(gff_index)
    def _read_gff_data(self, ref_id):
        gff_annotation = dict()
        for offset, size in self.gff_index.get(ref_id, list()):
            self.gff_data.seek(offset)
            for line in self.gff_data.read(size).strip("\n").split("\n"):
                if not line.startswith("#"):
                    line = line.strip().split("\t")
                    gff_annotation.setdefault((line[0], int(line[3]), int(line[4]) + 1), list()).append(line[8])
        if not gff_annotation:
            print("WARNING: contig {contig} does not have an annotation in the index.".format(contig=ref_id))
        return gff_annotation

    def process_cache(self):
        for qname, cached in self.read_cache.items():
            if len(cached) == 1:
                start, end = cached[0].start, cached[0].end
            elif len(cached) == 2:
                start, end = min(cached[0].start, cached[1].start), max(cached[0].end, cached[1].end)
            else:
                print("WARNING: more than two primary alignments for {qname} ({n}). Ignoring.".format(qname=qname, n=len(cached)))
                continue
            ovl_features = self.get_overlaps(start, end)
            print(*cached, sep="\n")
            print(*ovl_features, sep="\n")
        self.read_cache.clear()

    def collect_secondary_alignments(self, bamfile):
        t0 = time.time()
        bam = BamFile(bamfile)
        reads_with_secaln = set(
            aln.qname for aln in bam.get_alignments(
                required_flags=0x100, disallowed_flags=0x800
            )
        )
        bam = BamFile(bamfile) # need to implement rewind?
        sec_alignments = dict()
        for aln in bam.get_alignments(disallowed_flags=0x900):
            if aln.qname in reads_with_secaln:
                sec_alignments.setdefault(aln.qname, set()).add(
                    (aln.rid, aln.start, aln.end)
                )
        t1 = time.time()
        print("Collected information for {n_alignments} secondary alignments in {n_seconds:.3f}s.".format(
            n_alignments=len(sec_alignments), n_seconds=t1-t0)
        )
        missing = len(reads_with_secaln.difference(sec_alignments))
        if missing:
            print("{n_missing} secondary alignments don't have primary alignment in file.".format(n_missing=missing))

        return sec_alignments

    def process_bam(self, bamfile):
        sec_alignments = self.collect_secondary_alignments(bamfile)
        t0 = time.time()
        bam = BamFile(bamfile)
        self.current_ref = None
        self.gff_annotation = dict()
        self.read_cache = dict()
        sec_cache = dict()
        #only using enumerate here, so could break after n alignments for testing purposes
        for i, aln in enumerate(bam.get_alignments(disallowed_flags=0x800)):
            if not aln.is_primary():
                primaries = sec_alignments.get(aln.qname, set())
                if not primaries:
                    print("WARNING: could not find primary alignments for {aln_qname}".format(aln_qname=aln.qname))
                    continue
                aln_key = (aln.rid, aln.start, aln.end)
                if aln_key in primaries:
                    # there are sometimes secondary alignments that seem to be identical to the primary
                    continue
                if aln_key in sec_cache.get(aln.qname, set()):
                    # if r1 and r2 are identical and generate two separate individual alignments, disregard those
                    continue
                sec_cache.setdefault(aln.qname, set()).add(aln_key)

            elif aln.is_paired() and aln.rid == aln.rnext:
                # due to filtered mates (ngless!), need to keep track of read pairs to avoid dual counts
                # do not support paired information in secondary information currently
                self.read_cache.setdefault(aln.qname, list()).append(aln)

            ref = bam.get_reference(aln.rid)
            if ref != self.current_ref:
                # clear cache
                self.process_cache()
                sec_cache.clear()
                # load current reference
                self.current_ref = ref
                self.gff_annotation = self._read_gff_data(self.current_ref)
                intervals = sorted([key[1:] for key in self.gff_annotation])
                self.interval_tree = IntervalTree.from_tuples(intervals)

            print(aln)
            ovl_features = self.get_overlaps(aln.start, aln.end)
            print(*ovl_features, sep="\n")

        # clear cache
        self.process_cache()
        t1 = time.time()

        print("Processed {n_align} primary alignments in {n_seconds:.3f}s.".format(
            n_align=i, n_seconds=t1-t0)
        )


    def get_overlaps(self, start, end):
        overlaps = self.interval_tree[start:end]
        features = list()
        for ovl in overlaps:
            key = (self.current_ref, ovl.begin, ovl.end)
            features.append(self.gff_annotation[key])
        return features

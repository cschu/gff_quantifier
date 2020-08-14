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

    def process_bam(self, bamfile):
        bam = BamFile(bamfile)
        self.current_ref = None
        self.gff_annotation = dict()
        self.read_cache = dict()
        #only using enumerate here, so could break after n alignments for testing purposes
        for i, aln in enumerate(bam.get_alignments()):
            if not aln.is_primary():
                # TODO: deal with multimappers
                continue
            if aln.is_paired() and aln.rid == aln.rnext:
                self.read_cache.setdefault(aln.qname, list()).append(aln)
            # yield rid, self._references[rid][0], pos, len_seq
            ref = bam.get_reference(aln.rid)
            # yield BamAlignment(qname, flag, rid, pos, mapq, cigar, next_rid, next_pos, tlen, len_seq, tags)
            if ref != self.current_ref:
                # clear cache
                self.process_cache()
                # load current reference
                self.current_ref = ref
                self.gff_annotation = self._read_gff_data(self.current_ref)
                intervals = sorted([key[1:] for key in self.gff_annotation])
                self.interval_tree = IntervalTree.from_tuples(intervals)

            #overlaps = interval_tree[aln.start:aln.end]
            print(aln)
            ovl_features = self.get_overlaps(aln.start, aln.end)
            #print(overlaps)
            print(*ovl_features, sep="\n")

        # clear cache
        self.process_cache()

    def get_overlaps(self, start, end):
        overlaps = self.interval_tree[start:end]
        features = list()
        for ovl in overlaps:
            key = (self.current_ref, ovl.begin, ovl.end)
            features.append(self.gff_annotation[key])
        return features

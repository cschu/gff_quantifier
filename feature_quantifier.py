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

    def process_bam(self, bamfile):
        bam = BamFile(bamfile)
        current_ref = None
        gff_annotation = dict()
        #only using enumerate here, so could break after n alignments for testing purposes
        for i, aln in enumerate(bam.get_alignments()):
            #yield rid, self._references[rid][0], pos, len_seq
            if aln[1] != current_ref:
                current_ref = aln[1]
                gff_annotation = self._read_gff_data(current_ref)
                intervals = sorted([key[1:] for key in gff_annotation])
                # print(intervals)
                interval_tree = IntervalTree.from_tuples(intervals)
            overlaps = interval_tree[aln[2]:aln[2] + aln[3]]
            print(aln)
            print(overlaps)

            for ovl in overlaps:
                key = (current_ref, ovl.begin, ovl.end)
                print(gff_annotation[key])

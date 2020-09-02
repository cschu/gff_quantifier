import os
import time
import sys
import gzip
from collections import Counter
import json

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
    def _read_gff_data(self, ref_id, include_payload=False):
        gff_annotation = dict()
        for offset, size in self.gff_index.get(ref_id, list()):
            self.gff_data.seek(offset)
            for line in self.gff_data.read(size).strip("\n").split("\n"):
                if not line.startswith("#"):
                    line = line.strip().split("\t")
                    features = None
                    if include_payload:
                        features = dict((item.split("=")[0], item.split("=")[1].split(",")) for item in line[8].strip().split(";"))
                    key = (line[0], int(line[3]), int(line[4]) + 1)
                    gff_annotation[key] = features
                    #gff_annotation.setdefault((line[0], int(line[3]), int(line[4]) + 1), list()).append(line[8])
        if not gff_annotation:
            print("WARNING: contig {contig} does not have an annotation in the index.".format(contig=ref_id), file=sys.stderr, flush=True)
        return gff_annotation

    def process_cache(self):
        for qname, cached in self.read_cache.items():
            if len(cached) == 1:
                _, start, end = cached[0] #[0].start, cached[0].end
            elif len(cached) == 2:
                _, start1, end1 = cached[0]
                _, start2, end2 = cached[1]
                coords = sorted((start1, start2, end1, end2))
                start, end = coords[0], coords[-1]
            else:
                print("WARNING: more than two primary alignments for {qname} ({n}). Ignoring.".format(qname=qname, n=len(cached)), file=sys.stderr, flush=True)
                continue

            self.update_overlap_counter(self.primary_counter, start, end)

        self.read_cache.clear()

    def collect_secondary_alignments(self, bam):
        t0 = time.time()
        reads_with_secaln = set(
            aln.qname for aln in bam.get_alignments(
                required_flags=0x100, disallowed_flags=0x800
            )
        )
        bam.rewind()
        sec_alignments = dict()
        for i, aln in enumerate(bam.get_alignments(disallowed_flags=0x900), start=1):
            if aln.qname in reads_with_secaln:
                sec_alignments.setdefault(aln.qname, set()).add(
                    (aln.rid, aln.start, aln.end)
                )
        t1 = time.time()
        print("Collected information for {n_sec_alignments} secondary alignments ({size} bytes) in {n_seconds:.3f}s. (total: {n_alignments})".format(
            size=sys.getsizeof(sec_alignments), n_sec_alignments=len(sec_alignments), n_seconds=t1-t0, n_alignments=len(sec_alignments) + i), flush=True,
        )
        missing = len(reads_with_secaln.difference(sec_alignments))
        if missing:
            print("{n_missing} secondary alignments don't have primary alignment in file.".format(n_missing=missing), flush=True)

        return sec_alignments

    def process_bam(self, bamfile, out_prefix):
        bam = BamFile(bamfile)
        sec_alignments = self.collect_secondary_alignments(bam)
        bam.rewind()
        t0 = time.time()
        self.current_ref = None
        self.current_rid = None
        self.gff_annotation = dict()
        self.read_cache = dict()
        sec_cache = dict()
        self.primary_counter, self.secondary_counter = dict(), dict()

        #only using enumerate here, so could break after n alignments for testing purposes
        for i, aln in enumerate(bam.get_alignments(disallowed_flags=0x800)):
            if i % 10000 == 0:
                print("READ_CACHE SIZE IS {size} bytes (nreads={nreads})".format(size=sys.getsizeof(self.read_cache), nreads=len(self.read_cache)), flush=True)

            counter = self.primary_counter
            ref = bam.get_reference(aln.rid)
            start, end = aln.start, aln.end

            if not aln.is_primary():
                counter = self.secondary_counter
                primaries = sec_alignments.get(aln.qname, set())
                if not primaries:
                    print("WARNING: could not find primary alignments for {aln_qname}".format(aln_qname=aln.qname), flush=True, file=sys.stderr)
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
                mates = self.read_cache.setdefault(aln.qname, list())
                if mates:
                    mrid, mstart, mend = mates[0]
                    if aln.rnext != mrid:
                        print("WARNING: alignment {qname} seems to be corrupted: {aln1} {aln2}".format(qname=aln.qname, aln1=str(aln), aln2=str(mates[0])), flush=True, file=sys.stderr)
                        continue # i don't think this ever happens
                    else:
                        coords = sorted((aln.start, mstart, aln.end, mend))
                        start, end = coords[0], coords[-1]
                    self.read_cache.pop(aln.qname)
                else:
                    #self.read_cache.setdefault(aln.qname, list()).append((aln.rid, aln.start, aln.end))
                    mates.append((aln.rid, aln.start, aln.end))
                    continue # right?

            if ref != self.current_ref:
                self.current_rid = aln.rid
                if self.current_ref is not None:
                    # clear cache
                    self.process_cache()
                    sec_cache.clear()
                # load current reference
                self.current_ref = ref
                self.gff_annotation = self._read_gff_data(self.current_ref)
                intervals = sorted([key[1:] for key in self.gff_annotation])
                self.interval_tree = IntervalTree.from_tuples(intervals)

            self.update_overlap_counter(counter, start, end)

        # clear cache
        self.process_cache()
        t1 = time.time()

        print("Processed {n_align} primary alignments in {n_seconds:.3f}s.".format(
            n_align=i, n_seconds=t1-t0), flush=True
        )

        self.dump_overlap_counters(bam, out_prefix)


    def dump_overlap_counters(self, bam, out_prefix):
        print("Dumping overlap counters...", flush=True)
        with open(out_prefix + ".seqname.txt", "w") as out:
            for ref_id in set(self.primary_counter.get("seqname", Counter())).union(self.secondary_counter.get("seqname", Counter())):
                primary_count = self.primary_counter.get("seqname", Counter())[ref_id]
                secondary_count = self.secondary_counter.get("seqname", Counter())[ref_id]
                print(ref_id, bam.get_reference(ref_id), primary_count, secondary_count, flush=True, sep="\t", file=out)
        del self.primary_counter["seqname"]
        del self.secondary_counter["seqname"]
        print("Size primary counter: {size} bytes".format(size=sys.getsizeof(self.primary_counter)))
        print("Size secondary counter: {size} bytes".format(size=sys.getsizeof(self.secondary_counter)))
        with open(out_prefix + ".feature_counts.txt", "w") as out:
            dump_counter = [dict(), dict()]
            for ref_id in sorted(set(self.primary_counter).union(self.secondary_counter)):
                ref_name = bam.get_reference(ref_id)
                gff_annotation = self._read_gff_data(ref_name, include_payload=True)
                for (start, end) in set(self.primary_counter.get(ref_id, Counter())).union(self.secondary_counter.get(ref_id, Counter())):
                    primary_count = self.primary_counter.get(ref_id, Counter())[(start, end)]
                    secondary_count = self.secondary_counter.get(ref_id, Counter())[(start, end)]
                    for feature_type, values in gff_annotation.get((ref_name, start, end), dict()).items():
                        if feature_type != "ID":
                            dump_counter[0].setdefault(feature_type, Counter()).update({v: primary_count for v in values})
                            dump_counter[1].setdefault(feature_type, Counter()).update({v: secondary_count for v in values})
            json.dump(dump_counter, out, indent="\t")


    def update_overlap_counter(self, counter, start, end):
        counter.setdefault("seqname", Counter())[self.current_rid] += 1
        overlaps = self.interval_tree[start:end]
        for ovl in overlaps:
            counter.setdefault(self.current_rid, Counter())[(ovl.begin, ovl.end)] += 1

    def get_overlaps(self, start, end):
        overlaps = self.interval_tree[start:end]
        features = list()
        for ovl in overlaps:
            key = (self.current_ref, ovl.begin, ovl.end)
            features.append(self.gff_annotation[key])
        return features

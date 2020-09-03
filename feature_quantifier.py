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
        if not gff_annotation:
            print("WARNING: contig {contig} does not have an annotation in the index.".format(contig=ref_id), file=sys.stderr, flush=True)
        return gff_annotation

    def process_cache(self):
        for qname, cached in self.read_cache.items():
            if len(cached) == 1:
                _, start, end = cached[0]
            elif len(cached) == 2:
                _, start1, end1 = cached[0]
                _, start2, end2 = cached[1]
                start, end = FeatureQuantifier.calculate_fragment_borders(start1, end1, start2, end2)
            else:
                print("WARNING: more than two primary alignments for {qname} ({n}). Ignoring.".format(qname=qname, n=len(cached)), file=sys.stderr, flush=True)
                continue

            self.update_overlap_counter(self.unique_counter, start, end)

        self.read_cache.clear()

    def collect_multimappers(self, bam):
        t0 = time.time()
        reads_with_secaln = Counter(
            aln.qname for aln in bam.get_alignments(
                required_flags=0x100, disallowed_flags=0x800
            )
        )
        bam.rewind()
        multimappers = dict()
        for aln_count, aln in enumerate(bam.get_alignments(disallowed_flags=0x900), start=1):
            if aln.qname in reads_with_secaln:
                multimappers.setdefault(aln.qname, [reads_with_secaln[aln.qname], set()])[1].add(
                    (aln.rid, aln.start, aln.end)
                )
        t1 = time.time()
        print("Collected information for {n_multimappers} secondary alignments ({size} bytes) in {n_seconds:.3f}s. (total: {n_alignments})".format(
            size=sys.getsizeof(multimappers), n_multimappers=len(multimappers), n_seconds=t1-t0, n_alignments=len(multimappers) + aln_count), flush=True,
        )
        missing = len(set(reads_with_secaln).difference(multimappers))
        if missing:
            print("{n_missing} secondary alignments don't have primary alignment in file.".format(n_missing=missing), flush=True)

        return multimappers

    def process_bam(self, bamfile, out_prefix):
        bam = BamFile(bamfile)
        multimappers = self.collect_multimappers(bam)
        bam.rewind()
        t0 = time.time()
        self.current_ref = None
        self.current_rid = None
        self.gff_annotation = dict()
        self.read_cache = dict()
        multimap_cache = dict()
        self.unique_counter, self.multimap_counter = dict(), dict()

        for aln_count, aln in enumerate(bam.get_alignments(disallowed_flags=0x800)):

            counter = self.unique_counter
            ref = bam.get_reference(aln.rid)[0]
            start, end = aln.start, aln.end

            if ref != self.current_ref:
                self.current_rid = aln.rid
                if self.current_ref is not None:
                    # clear cache
                    self.process_cache()
                    # multimap_cache.clear()
                # load current reference
                self.current_ref = ref
                self.gff_annotation = self._read_gff_data(self.current_ref)
                intervals = sorted([key[1:] for key in self.gff_annotation])
                self.interval_tree = IntervalTree.from_tuples(intervals)

            mm_count, primaries = multimappers.get(aln.qname, [0, set()])
            if not aln.is_primary() or primaries:
                counter = self.multimap_counter
                if not primaries:
                    # ignore secondary alignments without existing primary alignment for now
                    print("WARNING: could not find primary alignments for {aln_qname}".format(aln_qname=aln.qname), flush=True, file=sys.stderr)
                    continue
                # we need to have the primary alignment flag in the key so we can process the primary alignment properly
                aln_key = (aln.rid, aln.start, aln.end, aln.is_primary())
                # special treatment for the secondary alignments
                # this is stuff i've seen coming out of ngless/bwa
                if not aln_key[-1]:
                    mm_count -= 1
                    if not mm_count and len(multimap_cache.get(aln.qname, set())) == len(primaries):
                        # all primary and secondary alignments of read have been seen
                        # TODO process alignments
                        multimap_cache.pop(aln.qname)
                        continue
                    multimappers[aln.qname][0] = mm_count

                    if aln_key[:-1] in primaries:
                        # there are sometimes secondary alignments that seem to be identical to the primary -> ignore
                        continue
                    if aln_key in multimap_cache.get(aln.qname, set()):
                        # if r1 and r2 are identical and generate two separate individual alignments -> ignore
                        # i don't know why that would happen: very short fragment sequenced twice, i.e. from each direction or r1/r2 mislabeled/duplicated?
                        continue
                    
                multimap_cache.setdefault(aln.qname, set()).add(aln_key)
                continue

            elif aln.is_paired() and aln.rid == aln.rnext:
                # need to keep track of read pairs to avoid dual counts
                # check if the mate has already been seen
                mates = self.read_cache.setdefault(aln.qname, list())
                if mates:
                    # if it has, calculate the total fragment size and remove the pair from the cache
                    if aln.rnext != mates[0][0]:
                        print("WARNING: alignment {qname} seems to be corrupted: {aln1} {aln2}".format(qname=aln.qname, aln1=str(aln), aln2=str(mates[0])), flush=True, file=sys.stderr)
                        continue # i don't think this ever happens
                    else:
                        start, end = FeatureQuantifier.calculate_fragment_borders(aln.start, aln.end, *mates[0][1:])
                    self.read_cache.pop(aln.qname)
                else:
                    # otherwise cache the first encountered mate and advance to the next read
                    mates.append((aln.rid, aln.start, aln.end))
                    continue

            # at this point only single-end reads and merged pairs should be processed here
            self.update_overlap_counter(counter, start, end)

        # clear cache
        self.process_cache()
        t1 = time.time()

        print("Processed {n_align} primary alignments in {n_seconds:.3f}s.".format(
            n_align=aln_count, n_seconds=t1-t0), flush=True
        )

        self.dump_overlap_counters(bam, out_prefix)
    
    @staticmethod
    def calculate_fragment_borders(start1, end1, start2, end2):
        coords = sorted((start1, end1, start2, end2))
        return coords[0], coords[-1]
    
        


    def dump_overlap_counters(self, bam, out_prefix):
        print("Dumping overlap counters...", flush=True)
        with open(out_prefix + ".seqname.txt", "w") as out:
            for ref_id in set(self.unique_counter.get("seqname", Counter())).union(self.multimap_counter.get("seqname", Counter())):
                primary_count = self.unique_counter.get("seqname", Counter())[ref_id]
                secondary_count = self.multimap_counter.get("seqname", Counter())[ref_id]
                print(ref_id, bam.get_reference(ref_id), primary_count, secondary_count, flush=True, sep="\t", file=out)
        del self.unique_counter["seqname"]
        try:
            del self.multimap_counter["seqname"]
        except:
            pass
        print("Size primary counter: {size} bytes".format(size=sys.getsizeof(self.unique_counter)))
        print("Size secondary counter: {size} bytes".format(size=sys.getsizeof(self.multimap_counter)))
        with open(out_prefix + ".feature_counts.txt", "w") as out:
            dump_counter = [dict(), dict()]
            for ref_id in sorted(set(self.unique_counter).union(self.multimap_counter)):
                ref_name = bam.get_reference(ref_id)[0]
                gff_annotation = self._read_gff_data(ref_name, include_payload=True)
                for (start, end) in set(self.unique_counter.get(ref_id, Counter())).union(self.multimap_counter.get(ref_id, Counter())):
                    primary_count = self.unique_counter.get(ref_id, Counter())[(start, end)]
                    secondary_count = self.multimap_counter.get(ref_id, Counter())[(start, end)]
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

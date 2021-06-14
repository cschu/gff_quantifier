import os
import time
import sys
from collections import Counter

import numpy


"""
normalizeCounts nmethod counts sizes
    | nmethod `elem` [NMScaled, NMFpkm] = do
        -- count vectors always include a -1 at this point (it is
        -- ignored in output if the user does not request it, but is
        -- always computed). Thus, we compute the sum without it and do
        -- not normalize it later:
        let totalCounts v = withVector v (VU.sum . VU.tail)
        initial <- totalCounts counts
        normalizeCounts NMNormed counts sizes
        afternorm <- totalCounts counts
        let factor
                | nmethod == NMScaled = initial / afternorm
                | otherwise = 1.0e9 / initial --- 1e6 [million fragments] * 1e3 [kilo basepairs] = 1e9
        liftIO $ forM_ [1.. VUM.length counts - 1] (VUM.unsafeModify counts (* factor))
"""

class OverlapCounter(dict):
    COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled"]

    @staticmethod
    def normalise_counts(counts, feature_len, scaling_factor):
        '''Returns raw, length-normalised, and scaled feature counts.'''
        normalised = counts / feature_len
        scaled = normalised * scaling_factor
        return counts, normalised, scaled
    def __init__(self, out_prefix, db, do_overlap_detection=True, strand_specific=False):
        self.out_prefix = out_prefix
        self.db = db
        self.seqcounts = Counter()
        self.ambig_seqcounts = Counter()
        self.featcounts = dict()
        self.unannotated_reads = 0
        self.ambig_counts = dict()
        self.has_ambig_counts = False
        self.feature_scaling_factors = dict()
        self.do_overlap_detection = do_overlap_detection
        self.strand_specific = strand_specific
        self.gene_counts = dict()
        self.coverage_intervals = dict()
        self.ambig_coverage = dict()

    def update_coverage_intervals(self, rid, intervals, coverage, ambig_aln=False):
        for (istart, iend), (cstart, cend) in zip(intervals, coverage):
            if ambig_aln:
                self.ambig_coverage.setdefault(rid, list()).append((rid, istart, iend, cstart, cend))
            else:
                self.coverage_intervals.setdefault(rid, dict()).setdefault((istart, iend), Counter())[(cstart, cend)] += 1

    def update_unique_counts(self, rid, aln=None, rev_strand=False):
        '''
        Generates a hit list from the overlaps resulting from an intervaltree query,
        adds the number of alternative alignments and stores the results for each reference sequence.
        '''
        if aln:
            overlaps, coverage = self.db.get_overlaps(*aln)
            if overlaps:
                self.setdefault(rid, Counter()).update((ovl.begin, ovl.end, rev_strand) for ovl in overlaps)
                self.update_coverage_intervals(rid, overlaps, coverage)
            else:
                self.unannotated_reads += 1
        if self.strand_specific and not self.do_overlap_detection:
            self.seqcounts[(rid, rev_strand)] += 1
        else:
            self.seqcounts[rid] += 1

    def _compute_count_vector(self, bins, rid, aln, strand, strand_specific):
        # calculate the count vector for the current region
        # count vectors are of the format (raw_uniq, normed_uniq, raw_ambi, normed_ambi)
        # att: _ambi counts are unique + ambiguous
        counts = numpy.zeros(bins)
        region_length = aln[1] - aln[0] + 1
        counts[0] = counts[1] = self.get(rid, dict()).get(aln, 0.0)
        counts[1] /= region_length
        counts[2] = counts[3] = counts[0] + self.ambig_counts.get(rid, dict()).get(aln, 0.0)
        counts[3] /= region_length
        if strand_specific:
            rev_strand = aln[-1]
            is_antisense = (strand == "+" and rev_strand) or (strand == "-" and not rev_strand)
            bin_offset = 8 if is_antisense else 4
            counts[bin_offset:bin_offset + 4] += counts[:4]

        return counts

    def _iterate_database(self, bins, bam, strand_specific):
        total_counts, feature_count_sums = numpy.zeros(4), dict()

        for ref, region_annotation in self.db.iterate():
            rid = bam.revlookup_reference(ref)
            if rid is not None:
                _, region_length = bam.get_reference(rid)
                counts = numpy.zeros(bins)

                if strand_specific:
                    counts[0] = counts[1] = self.seqcounts[(rid, True)] + self.seqcounts[(rid, False)]
                    counts[1] /= region_length
                    counts[2] = counts[3] = counts[0] + self.ambig_seqcounts[(rid, True)] + self.ambig_seqcounts[(rid, False)]
                    counts[3] /= region_length
                    counts[4] = counts[5] = self.seqcounts[(rid, True)]
                    counts[5] /= region_length
                    counts[6] = counts[7] = counts[4] + self.ambig_seqcounts[(rid, True)]
                    counts[7] /= region_length
                    counts[8] = counts[9] = self.seqcounts[(rid, False)]
                    counts[9] /= region_length
                    counts[10] = counts[11] = counts[8] + self.ambig_seqcounts[(rid, False)]
                    counts[11] /= region_length
                else:
                    counts[0] = counts[1] = self.seqcounts[rid]
                    counts[1] /= region_length
                    counts[2] = counts[3] = counts[0] + self.ambig_seqcounts[rid]
                    counts[3] /= region_length

                total_counts, feature_count_sums = self._distribute_feature_counts(bins, counts, region_annotation[1:], total_counts, feature_count_sums)
                self._add_count_vector(counts, ref, self.gene_counts, bins)

        return total_counts, feature_count_sums


    def _iterate_bedcounts(self, bins, feature_lengths, strand_specific):
        total_counts, feature_count_sums = numpy.zeros(4), dict()

        for rid in set(self.seqcounts).union(self.ambig_seqcounts):
            region_length = feature_lengths.get(rid)
            feature = rid.split("::")[-1]

            if region_length is None:
                raise ValueError(f"Cannot determine length of reference {rid}.")
            counts = numpy.zeros(bins)

            if strand_specific:
                counts[0] = counts[1] = self.seqcounts[(rid, True)] + self.seqcounts[(rid, False)]
                counts[1] /= region_length
                counts[2] = counts[3] = counts[0] + self.ambig_seqcounts[(rid, True)] + self.ambig_seqcounts[(rid, False)]
                counts[3] /= region_length
                counts[4] = counts[5] = self.seqcounts[(rid, True)]
                counts[5] /= region_length
                counts[6] = counts[7] = counts[4] + self.ambig_seqcounts[(rid, True)]
                counts[7] /= region_length
                counts[8] = counts[9] = self.seqcounts[(rid, False)]
                counts[9] /= region_length
                counts[10] = counts[11] = counts[8] + self.ambig_seqcounts[(rid, False)]
                counts[11] /= region_length
            else:
                counts[0] = counts[1] = self.seqcounts[rid]
                counts[1] /= region_length
                counts[2] = counts[3] = counts[0] + self.ambig_seqcounts[rid]
                counts[3] /= region_length

            total_counts, feature_count_sums = self._distribute_feature_counts(bins, counts, [("feature", (feature,))], total_counts, feature_count_sums)

        return total_counts, feature_count_sums


    def _distribute_feature_counts(self, bins, counts, region_annotation, total_counts, feature_count_sums):
        for ftype, ftype_counts in region_annotation:
            total_fcounts = feature_count_sums.setdefault(ftype, numpy.zeros(4))
            for i, ft_ct in enumerate(ftype_counts, start=1):
                fcounts = self.featcounts.setdefault(ftype, dict())
                self._add_count_vector(counts, ft_ct, fcounts, bins)

            inc = counts[:4] * i
            total_counts += inc
            total_fcounts += inc

        return total_counts, feature_count_sums

    def _add_count_vector(self, count_vector, target_id, target_ctr, bins):
        counts = target_ctr.setdefault(target_id, numpy.zeros(bins))
        counts += count_vector

    def _iterate_counts(self, bins, bam, strand_specific):
        total_counts, feature_count_sums = numpy.zeros(4), dict()

        for rid in set(self.keys()).union(self.ambig_counts):
            ref = bam.get_reference(rid)[0]
            for start, end, rev_strand in set(self.get(rid, set())).union(self.ambig_counts.get(rid, set())):
                # region_annotation is a tuple of key-value pairs: (strand, func_category1: subcategories, func_category2: subcategories, ...)
                # the first is the strand, the second is the gene id, the rest are the features
                region_annotation = self.db.get_data(ref, start, end)
                counts = self._compute_count_vector(bins, rid, (start, end, rev_strand), region_annotation[0][1], strand_specific)
                # how to extract the gene counts in genome mode?
                self._add_count_vector(counts, region_annotation[1][1], self.gene_counts, bins)

                # distribute the counts to the associated functional (sub-)categories
                total_counts, feature_count_sums = self._distribute_feature_counts(bins, counts, region_annotation[2:], total_counts, feature_count_sums)

        return total_counts, feature_count_sums


    def annotate_counts(self, bamfile=None, feature_lengths=None, itermode="counts"):
        """
        Distributes read counts against aligned regions to the associated functional categories.
        """
        print("Processing counts ...", flush=True)
        t0 = time.time()
        bins = 12 if self.strand_specific else 4


        annotation_f = {
            "counts": self._iterate_counts,
            "database": self._iterate_database,
            "bedcounts": self._iterate_bedcounts
        }.get(itermode)

        if not annotation_f:
            raise ValueError(f"Unknown annotation_mode {itermode}")

        if bool(bamfile) == bool(feature_lengths):
            raise ValueError(f"Cannot determine feature lengths.")
        if not feature_lengths:
            feature_lengths = bamfile

        total_counts, feature_count_sums = annotation_f(bins, feature_lengths, self.strand_specific)

        # calculate the scaling factors
        total, total_normed, total_ambi, total_ambi_normed = total_counts

        default_scaling_factor = 0
        self.feature_scaling_factors["total"] = (total / total_normed) if total_normed else default_scaling_factor
        self.feature_scaling_factors["total_ambi"] = (total_ambi / total_ambi_normed) if total_ambi_normed else default_scaling_factor

        for ftype, counts in feature_count_sums.items():
            total, total_normed, total_ambi, total_ambi_normed = counts
            self.feature_scaling_factors[ftype] = (
                (total / total_normed) if total_normed else default_scaling_factor,
                (total_ambi / total_ambi_normed) if total_ambi_normed else default_scaling_factor
            )

        self.clear()
        self.ambig_counts.clear()
        t1 = time.time()
        print("Processed counts in {n_seconds}s.".format(n_seconds=t1-t0), flush=True)


    def update_ambiguous_counts(self, hits, n_aln, unannotated=0, bam=None, feat_distmode="all1"):
        self.has_ambig_counts = True
        self.unannotated_reads += unannotated
        if self.strand_specific and not self.do_overlap_detection:
            n_total = sum(self.seqcounts[(rid, True)] + self.seqcounts[(rid, False)] for rid in hits)
        else:
            n_total = sum(self.seqcounts[rid] for rid in hits)
        for rid, regions in hits.items():
            if self.do_overlap_detection:
                seen_regions = set()
                for start, end, cstart, cend, rev_str in regions:
                    if region[:3] not in seen_regions:
                        self.ambig_counts.setdefault(rid, Counter())[(start, end, rev_str)] += (1 / n_aln) if feat_distmode == "1overN" else 1

                    self.update_coverage_intervals(rid, [(start, end)], [(cstart, cend)], ambig_aln=True)
                    seen_regions.add(region[:3])

                if n_total and self.seqcounts[rid]:
                    self.ambig_seqcounts[rid] += self.seqcounts[rid] / n_total * len(hits)
                else:
                    self.ambig_seqcounts[rid] += 1 / len(hits)
            else:
                start, end, rev_str, cstart, cend = list(regions)[0]
                key = (rid, rev_str) if self.strand_specific else rid
                self.ambig_seqcounts[key] += (1 / n_aln) if feat_distmode == "1overN" else 1
                self.update_coverage_intervals(rid, [(start, end)], [(cstart, cend)], ambig_aln=True)

    @staticmethod
    def calculate_seqcount_scaling_factor(counts, bam):
        raw_total = sum(counts.values())
        normed_total = sum(count / bam.get_reference(rid)[1] for rid, count in counts.items())
        return raw_total / normed_total
    @staticmethod
    def calculate_feature_scaling_factor(counts, include_ambig=False):
        raw_total, normed_total = 0, 0
        for raw, norm, raw_ambi, norm_ambi, *_ in counts.values():
            raw_total += raw
            normed_total += norm
            if include_ambig:
                raw_total += raw_ambi
                normed_total += norm_ambi
        return raw_total / normed_total


    def _compile_output_row(self, counts, scaling_factor=1, ambig_scaling_factor=1):
        p, row = 0, list()
        # unique counts
        row.extend(counts[p:p + 2])
        row.append(row[-1] * scaling_factor)
        p += 2
        # ambiguous counts
        if self.has_ambig_counts:
            row.extend(counts[p:p + 2])
            row.append(row[-1] * ambig_scaling_factor)
            p += 2
        # sense-strand unique
        if self.strand_specific:
            row.extend(counts[p:p + 2])
            row.append(row[-1] * scaling_factor)
            p += 2
            # sense-strand ambiguous
            if self.has_ambig_counts:
                row.extend(counts[p:p + 2])
                row.append(row[-1] * ambig_scaling_factor)
                p += 2
            # antisense-strand unique
            row.extend(counts[p:p + 2])
            row.append(row[-1] * scaling_factor)
            p += 2
            if self.has_ambig_counts:
                row.extend(counts[p:p + 2])
                row.append(row[-1] * ambig_scaling_factor)
        return row

    def get_header(self):
        header = list()
        header.extend(f"uniq_{element}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS)
        if self.has_ambig_counts:
            header.extend(f"combined_{element}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS)
        if self.strand_specific:
            for strand in ("ss", "as"):
                header.extend(f"uniq_{element}_{strand}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS)
                if self.has_ambig_counts:
                    header.extend(f"combined_{element}_{strand}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS)
        return header


    def dump_counts(self, bam=None):

        SEQ_COUNT_HEADER = ["seqid_int", "seqid", "length"] + OverlapCounter.COUNT_HEADER_ELEMENTS
        counts_template = "{:.5f}\t{:.5f}\t{:.5f}"
        header = self.get_header()

        print("Dumping overlap counters...", flush=True)
        print("Has ambiguous counts:", self.has_ambig_counts, flush=True)
        if self.strand_specific and not self.do_overlap_detection:
             _seqcounts = Counter()
             for (rid, rev_strand), count in self.seqcounts.items():
                 _seqcounts[rid] += count
             self.seqcounts = _seqcounts
             if self.has_ambig_counts:
                 _seqcounts = Counter()
                 for (rid, rev_strand), count in self.ambig_seqcounts.items():
                     _seqcounts[rid] += count
                 self.ambig_seqcounts = _seqcounts

        if self.gene_counts:
            with open(f"{self.out_prefix}.gene_counts.txt", "w") as gene_out:
                print("gene", *header, sep="\t", file=gene_out, flush=True)
                for gene, g_counts in sorted(self.gene_counts.items()):
                    out_row = self._compile_output_row(
                        g_counts,
                        scaling_factor=self.feature_scaling_factors["total"],
                        ambig_scaling_factor=self.feature_scaling_factors["total_ambi"]
                    )
                    if isinstance(gene, tuple):
+                       gene = gene[0]
                    print(gene, out_row[0], *(f"{c:.5f}" for c in out_row[1:]), flush=True, sep="\t", file=gene_out)

        if bam:
            with open(f"{self.out_prefix}.seqname.uniq.txt", "w") as seq_out:
                print(*SEQ_COUNT_HEADER, sep="\t", flush=True, file=seq_out)
                if sum(self.seqcounts.values()):
                    seqcount_scaling_factor = OverlapCounter.calculate_seqcount_scaling_factor(self.seqcounts, bam)
                    for rid, count in self.seqcounts.items():
                        seq_id, seq_len = bam.get_reference(rid)
                        print(rid, seq_id, seq_len, count, "{:.5f}\t{:.5f}".format(*OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)[1:]), flush=True, sep="\t", file=seq_out)

        with open(f"{self.out_prefix}.feature_counts.txt", "w") as feat_out:
            print("subfeature", *header, sep="\t", file=feat_out, flush=True)
            print("unannotated", self.unannotated_reads, sep="\t", file=feat_out, flush=True)
            for ftype, counts in sorted(self.featcounts.items()):
                print(f"#{ftype}", file=feat_out, flush=True)
                scaling_factor, ambig_scaling_factor = self.feature_scaling_factors[ftype]
                for subf, sf_counts in sorted(counts.items()):
                    out_row = self._compile_output_row(sf_counts, scaling_factor=scaling_factor, ambig_scaling_factor=ambig_scaling_factor)
                    print(subf, out_row[0], *(f"{c:.5f}" for c in out_row[1:]), flush=True, sep="\t", file=feat_out)

        if bam and self.ambig_seqcounts:
                with open(f"{self.out_prefix}.seqname.dist1.txt", "w") as seq_out:
                    print(*SEQ_COUNT_HEADER, sep="\t", flush=True, file=seq_out)
                    self.seqcounts.update(self.ambig_seqcounts)
                    seqcount_scaling_factor = OverlapCounter.calculate_seqcount_scaling_factor(self.seqcounts, bam)
                    for rid, count in self.seqcounts.items():
                        seq_id, seq_len = bam.get_reference(rid)
                        print(rid, seq_id, seq_len, counts_template.format(*OverlapCounter.normalise_counts(count, seq_len, seqcount_scaling_factor)), flush=True, sep="\t", file=seq_out)

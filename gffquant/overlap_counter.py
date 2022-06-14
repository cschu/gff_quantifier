# pylint: skip-file

import gzip
import time
from collections import Counter

import numpy as np


"""
WARNING: FILE IS OBSOLETE, ONLY KEPT FOR RESTORING COVERAGE-FUNCTIONALITY
"""


class OverlapCounter(dict):
    COUNT_HEADER_ELEMENTS = ["raw", "lnorm", "scaled"]

    @staticmethod
    def normalise_counts(counts, feature_len, scaling_factor):
        """Returns raw, length-normalised, and scaled feature counts."""
        normalised = counts / feature_len
        scaled = normalised * scaling_factor
        return counts, normalised, scaled

    def __init__(
        self,
        out_prefix,
        db,
        do_overlap_detection=True,
        strand_specific=False,
        feature_distribution="1overN",
    ):
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
        self.feature_distribution = feature_distribution

    def update_unique_counts(self, count_stream):
        strand_specific = self.strand_specific and not self.do_overlap_detection
        for counts, aln_count, unaligned in count_stream:
            for rid, hits in counts.items():
                strands = set()
                for ostart, oend, rev_strand, cstart, cend in hits:
                    if ostart is not None:
                        self.setdefault(rid, Counter())[(ostart, oend, rev_strand)] += 1

                        if self.db.reference_type == "domain":
                            self.coverage_intervals.setdefault(rid, dict()).setdefault(
                                (ostart, oend), Counter()
                            )[(cstart, cend)] += 1
                    strands.add(rev_strand)

                if strand_specific:
                    for strand in strands:
                        seqcount_key = (rid, rev_strand)
                        self.seqcounts[
                            seqcount_key
                        ] += 1  # aln_count  # this overcounts reads when they overlap multiple features
                else:
                    self.seqcounts[rid] += 1

            self.unannotated_reads += unaligned

    def update_ambiguous_counts(self, count_stream):
        count_stream = tuple(
            count_stream
        )  # we need to access stream-length before processing individual items. do we, though?
        self.has_ambig_counts = True
        strand_specific = self.strand_specific and not self.do_overlap_detection

        for counts, aln_count, unaligned in count_stream:
            coverage = {}
            self.unannotated_reads += unaligned
            increment = (
                (1 / aln_count) if self.feature_distribution == "1overN" else 1
            )  # 1overN = lavern. Maya <3
            if strand_specific:
                n_total = sum(
                    self.seqcounts[(rid, True)] + self.seqcounts[(rid, False)]
                    for rid in counts
                )
            else:
                n_total = sum(self.seqcounts[rid] for rid in counts)

            # gene mode counts come as                 hits = {(None, None, rev_strand, None, None)}
            # ovl mode counts come as hits = {(ostart, oend, rev_strand, cstart, cend)}
            for rid, hits in counts.items():
                regions, strands = {}, set()
                for ostart, oend, rev_strand, cstart, cend in hits:
                    if self.do_overlap_detection:  # ostart is None or ostart == -1:
                        regions.setdefault((ostart, oend, rev_strand), set()).add(
                            (cstart, cend)
                        )
                    strands.add(rev_strand)

                # pull in the overlap data
                for region, covered in regions.items():
                    self.ambig_counts.setdefault(rid, Counter())[region] += increment

                    ostart, oend, rev_strand = region

                    # coverage requires an overlap-mode
                    # (coordinates are not collected in gene mode), memory issue!
                    if self.db.reference_type == "domain":
                        for cstart, cend in covered:
                            coverage.setdefault(rid, dict()).setdefault(
                                (ostart, oend), list()
                            ).append((cstart, cend))

                if strand_specific:
                    for strand in strands:
                        seqcount_key = (rid, strand)
                        if n_total and self.seqcounts[seqcount_key]:
                            seq_increment = (
                                self.seqcounts[seqcount_key] / n_total * aln_count
                            )
                        else:
                            seq_increment = 1 / aln_count
                        self.ambig_seqcounts[seqcount_key] += seq_increment

                else:
                    seqcount_key = rid
                    if n_total and self.seqcounts[seqcount_key]:
                        seq_increment = (
                            self.seqcounts[seqcount_key] / n_total * aln_count
                        )
                    else:
                        seq_increment = 1 / aln_count
                    self.ambig_seqcounts[seqcount_key] += seq_increment

                """
				# seqcounts via strands
				# in gene mode feature counts are derived from seqcounts
				for strand in strands:
					seqcount_key = (rid, strand) if self.strand_specific else rid
					if n_total and self.seqcounts[seqcount_key]:
						seq_increment = self.seqcounts[seqcount_key] / n_total * aln_count
					else:
						seq_increment = 1 / aln_count

					self.ambig_seqcounts[seqcount_key] += seq_increment
				"""

            self.update_ambig_coverage(coverage, aln_count)

    def update_ambig_coverage(self, coverage_data, n_aln):
        pseudocount = 1e-10
        # this is data from one(!) alignment group
        # coverage_data.setdefault(rid, dict()).setdefault((start, end), list()).append((cstart, cend))
        rev_strand = False
        dist1_coverage = Counter()
        for rid, regions in coverage_data.items():
            for (start, end, *_), overlaps in regions.items():
                dist1_coverage[(rid, start, end)] = self.get(rid, Counter()).get(
                    (start, end, rev_strand), 0.0
                )
        n_uniq = sum(dist1_coverage.values())
        for rid, regions in coverage_data.items():
            for (start, end, *_), overlaps in regions.items():
                for cstart, cend in overlaps:
                    cov_interval = self.ambig_coverage.setdefault(
                        rid, dict()
                    ).setdefault((start, end, rev_strand), Counter())
                    increment = (
                        (pseudocount + dist1_coverage[(rid, start, end)] / n_uniq)
                        if n_uniq
                        else (pseudocount + 1 / n_aln)
                    )
                    cov_interval[(cstart, cend)] += increment

    def _compute_count_vector(self, bins, rid, aln, strand, strand_specific):
        # calculate the count vector for the current region
        # count vectors are of the format (raw_uniq, normed_uniq, raw_ambi, normed_ambi)
        # att: _ambi counts are unique + ambiguous
        counts = np.zeros(bins)
        region_length = aln[1] - aln[0] + 1
        counts[0] = counts[1] = self.get(rid, dict()).get(aln, 0.0)
        counts[1] /= region_length
        counts[2] = counts[3] = counts[0] + self.ambig_counts.get(rid, dict()).get(
            aln, 0.0
        )
        counts[3] /= region_length
        if strand_specific:
            rev_strand = aln[-1]
            is_antisense = (strand == "+" and rev_strand) or (
                strand == "-" and not rev_strand
            )
            bin_offset = 8 if is_antisense else 4
            counts[bin_offset:bin_offset + 4] += counts[:4]

        return counts

    def _compute_genes_count_vector(self, rid, length, strand_specific=False):
        counts = np.zeros(12 if strand_specific else 4)
        if strand_specific:
            PLUS_STRAND, MINUS_STRAND = True, False
            uniq_plus = self.seqcounts[(rid, PLUS_STRAND)]
            uniq_minus = self.seqcounts[(rid, MINUS_STRAND)]
            ambig_plus = self.ambig_seqcounts[(rid, PLUS_STRAND)]
            ambig_minus = self.ambig_seqcounts[(rid, MINUS_STRAND)]
            counts[0] = counts[1] = uniq_plus + uniq_minus
            counts[1] /= length
            counts[2] = counts[3] = counts[0] + ambig_plus + ambig_minus
            counts[3] /= length
            counts[4] = counts[5] = uniq_plus
            counts[5] /= length
            counts[6] = counts[7] = counts[4] + ambig_plus
            counts[7] /= length
            counts[8] = counts[9] = uniq_minus
            counts[9] /= length
            counts[10] = counts[11] = counts[8] + ambig_minus
            counts[11] /= length
        else:
            counts[0] = counts[1] = self.seqcounts[rid]
            counts[1] /= length
            counts[2] = counts[3] = counts[0] + self.ambig_seqcounts[rid]
            counts[3] /= length
        return counts

    def _iterate_database(self, bins, bam, strand_specific):
        total_counts, feature_count_sums = np.zeros(4), dict()

        for ref, region_annotation in self.db.iterate():
            rid = bam.revlookup_reference(ref)
            if rid is not None:
                _, region_length = bam.get_reference(rid)
                counts = self._compute_genes_count_vector(
                    rid, region_length, strand_specific=strand_specific
                )
                total_counts, feature_count_sums = self._distribute_feature_counts(
                    bins,
                    counts,
                    region_annotation[1:],
                    total_counts,
                    feature_count_sums,
                )
                self._add_count_vector(counts, ref, self.gene_counts, bins)

        return total_counts, feature_count_sums

    def _distribute_feature_counts(
        self, bins, counts, region_annotation, total_counts, feature_count_sums
    ):
        for ftype, ftype_counts in region_annotation:
            total_fcounts = feature_count_sums.setdefault(ftype, np.zeros(4))
            for i, ft_ct in enumerate(ftype_counts, start=1):
                fcounts = self.featcounts.setdefault(ftype, dict())
                self._add_count_vector(counts, ft_ct, fcounts, bins)

            try:
                inc = counts[:4] * i
            except NameError:
                inc = 0

            total_counts += inc
            total_fcounts += inc

        return total_counts, feature_count_sums

    def _add_count_vector(self, count_vector, target_id, target_ctr, bins):
        counts = target_ctr.setdefault(target_id, np.zeros(bins))
        counts += count_vector

    def _iterate_counts(self, bins, bam, strand_specific):
        total_counts, feature_count_sums = np.zeros(4), dict()

        for rid in set(self.keys()).union(self.ambig_counts):
            ref = bam.get_reference(rid[0] if isinstance(rid, tuple) else rid)[0]
            regions = set(self.get(rid, set())).union(self.ambig_counts.get(rid, set()))
            for start, end, rev_strand in regions:
                # region_annotation is a tuple of key-value pairs:
                # (strand, func_category1: subcategories, func_category2: subcategories, ...)
                # the first is the strand, the second is the gene id, the rest are the features
                region_annotation = self.db.get_data(ref, start, end)
                counts = self._compute_count_vector(
                    bins,
                    rid,
                    (start, end, rev_strand),
                    region_annotation[0][1],
                    strand_specific,
                )
                # how to extract the gene counts in genome mode?
                self._add_count_vector(
                    counts, region_annotation[1][1], self.gene_counts, bins
                )

                # distribute the counts to the associated functional (sub-)categories
                total_counts, feature_count_sums = self._distribute_feature_counts(
                    bins,
                    counts,
                    region_annotation[2:],
                    total_counts,
                    feature_count_sums,
                )

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
        }.get(itermode)

        if not annotation_f:
            raise ValueError(f"Unknown annotation_mode {itermode}")

        if bool(bamfile) == bool(feature_lengths):
            raise ValueError("Cannot determine feature lengths.")
        if not feature_lengths:
            feature_lengths = bamfile

        total_counts, feature_count_sums = annotation_f(
            bins, feature_lengths, self.strand_specific
        )

        # calculate the scaling factors
        total, total_normed, total_ambi, total_ambi_normed = total_counts

        default_scaling_factor = 0
        self.feature_scaling_factors = {
            "total": (total / total_normed) if total_normed else default_scaling_factor,
            "total_ambi": (total_ambi / total_ambi_normed)
            if total_ambi_normed
            else default_scaling_factor,
        }

        for ftype, counts in feature_count_sums.items():
            total, total_normed, total_ambi, total_ambi_normed = counts
            self.feature_scaling_factors[ftype] = (
                (total / total_normed) if total_normed else default_scaling_factor,
                (total_ambi / total_ambi_normed)
                if total_ambi_normed
                else default_scaling_factor,
            )

        self.clear()
        self.ambig_counts.clear()
        t1 = time.time()
        print("Processed counts in {n_seconds}s.".format(n_seconds=t1 - t0), flush=True)

    @staticmethod
    def calculate_seqcount_scaling_factor(counts, bam):
        raw_total = sum(counts.values())
        normed_total = sum(
            count / bam.get_reference(rid[0] if isinstance(rid, tuple) else rid)[1]
            for rid, count in counts.items()
        )
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
        # unique counts
        row.extend(counts[p:p + 2])
        row.append(row[-1] * scaling_factor)
        p += 2
        # ambiguous counts
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
        header += (
            f"uniq_{element}" for element in OverlapCounter.COUNT_HEADER_ELEMENTS
        )
        if self.has_ambig_counts:
            header += (
                f"combined_{element}"
                for element in OverlapCounter.COUNT_HEADER_ELEMENTS
            )
        if self.strand_specific:
            for strand in ("ss", "as"):
                header += (
                    f"uniq_{element}_{strand}"
                    for element in OverlapCounter.COUNT_HEADER_ELEMENTS
                )
                if self.has_ambig_counts:
                    header += (
                        f"combined_{element}_{strand}"
                        for element in OverlapCounter.COUNT_HEADER_ELEMENTS
                    )
        return header

    def _dump_feature_counts(self):
        for ftype, counts in sorted(self.featcounts.items()):
            with gzip.open(f"{self.out_prefix}.{ftype}.txt.gz", "wt") as feat_out:
                print("feature", *self.get_header(), sep="\t", file=feat_out, flush=True)
                print("unannotated", self.unannotated_reads, sep="\t", file=feat_out, flush=True)
                scaling_factor, ambig_scaling_factor = self.feature_scaling_factors[ftype]
                for subf, sf_counts in sorted(counts.items()):
                    out_row = self._compile_output_row(
                        sf_counts, scaling_factor=scaling_factor, ambig_scaling_factor=ambig_scaling_factor
                    )
                    print(
                        subf, *(f"{c:.5f}" for c in out_row),
                        flush=True, sep="\t", file=feat_out
                    )

    def _dump_seq_counts(self, bam):
        SEQ_COUNT_HEADER = [
            "seqid_int",
            "seqid",
            "length",
        ] + OverlapCounter.COUNT_HEADER_ELEMENTS

        if self.strand_specific and not self.do_overlap_detection:
            _seqcounts = Counter()
            for (rid, _), count in self.seqcounts.items():
                _seqcounts[rid] += count
            self.seqcounts = _seqcounts
            if self.has_ambig_counts:
                _seqcounts = Counter()
                for (rid, _), count in self.ambig_seqcounts.items():
                    _seqcounts[rid] += count
                self.ambig_seqcounts = _seqcounts

        with gzip.open(f"{self.out_prefix}.seqname.uniq.txt.gz", "wt") as seq_out:
            print(*SEQ_COUNT_HEADER, sep="\t", flush=True, file=seq_out)
            if sum(self.seqcounts.values()):
                seqcount_scaling_factor = (
                    OverlapCounter.calculate_seqcount_scaling_factor(
                        self.seqcounts, bam
                    )
                )
                for rid, count in self.seqcounts.items():
                    seq_id, seq_len = bam.get_reference(
                        rid[0] if isinstance(rid, tuple) else rid
                    )
                    norm, scaled = OverlapCounter.normalise_counts(
                        count, seq_len, seqcount_scaling_factor
                    )[1:]
                    print(
                        rid,
                        seq_id,
                        seq_len,
                        count,
                        f"{norm:.5f}",
                        f"{scaled:.5f}",
                        flush=True,
                        sep="\t",
                        file=seq_out,
                    )

        if self.ambig_seqcounts:
            with gzip.open(f"{self.out_prefix}.seqname.dist1.txt.gz", "wt") as seq_out:
                print(*SEQ_COUNT_HEADER, sep="\t", flush=True, file=seq_out)
                self.seqcounts.update(self.ambig_seqcounts)
                seqcount_scaling_factor = (
                    OverlapCounter.calculate_seqcount_scaling_factor(
                        self.seqcounts, bam
                    )
                )
                for rid, count in self.seqcounts.items():
                    seq_id, seq_len = bam.get_reference(
                        rid[0] if isinstance(rid, tuple) else rid
                    )
                    raw, norm, scaled = OverlapCounter.normalise_counts(
                        count, seq_len, seqcount_scaling_factor
                    )
                    print(
                        rid,
                        seq_id,
                        seq_len,
                        f"{raw:.5f}",
                        f"{norm:.5f}",
                        f"{scaled:.5f}",
                        flush=True,
                        sep="\t",
                        file=seq_out,
                    )

    def _dump_gene_counts(self, bam=None):
        with gzip.open(f"{self.out_prefix}.gene_counts.txt.gz", "wt") as gene_out:
            print("gene", *self.get_header(), sep="\t", file=gene_out, flush=True)
            if self.do_overlap_detection:
                gene_counts = self.gene_counts
                scaling_factor = self.feature_scaling_factors["total"]
                ambig_scaling_factor = self.feature_scaling_factors["total_ambi"]
            else:
                gene_ids = set(
                    (gene_id[0] if isinstance(gene_id, tuple) else gene_id)
                    for gene_id in set(self.seqcounts).union(self.ambig_seqcounts)
                )
                gene_counts = {
                    gene_id: self._compute_genes_count_vector(
                        gene_id,
                        bam.get_reference(gene_id[0] if isinstance(gene_id, tuple) else gene_id)[1],
                        strand_specific=self.strand_specific,
                    )
                    for gene_id in gene_ids
                }
                counts_for_scaling = np.zeros(4)
                for counts in gene_counts.values():
                    counts_for_scaling += counts[:4]
                scaling_factor = counts_for_scaling[0] / counts_for_scaling[1]
                ambig_scaling_factor = counts_for_scaling[2] / counts_for_scaling[3]

            for gene, g_counts in sorted(gene_counts.items()):
                out_row = self._compile_output_row(
                    g_counts,
                    scaling_factor=scaling_factor,
                    ambig_scaling_factor=ambig_scaling_factor,
                )
                if not self.do_overlap_detection:
                    gene, _ = bam.get_reference(gene[0] if isinstance(gene, tuple) else gene)
                print(
                    gene,
                    out_row[0],
                    *(f"{c:.5f}" for c in out_row[1:]),
                    flush=True,
                    sep="\t",
                    file=gene_out,
                )

    def summarise_coverage(self, bam):
        ref_coverage, ambig_ref_coverage = dict(), dict()
        pos_feature_dict = dict()

        for rid, intervals in self.coverage_intervals.items():
            ref, _ = bam.get_reference(rid[0] if isinstance(rid, tuple) else rid)
            print("REF", ref, rid)
            for (start, end), overlaps in intervals.items():
                print(start, end, overlaps)
                features = self.db.db.get(ref, dict()).get((start, end), list())
                print("FEATURES", features)
                region = (rid, start, end)
                print("REGION", region)
                print("OVERLAPS", overlaps)
                region_coverage = Counter({p: 0 for p in range(start, end + 1)})
                ref_coverage.setdefault(region, Counter()).update(region_coverage)
                print("refcov", region, ref_coverage[region])
                ambig_ref_coverage.setdefault(region, Counter()).update(region_coverage)
                pos_feature_dict.setdefault(region, set()).update(features)
                for (ovl_start, ovl_end), count in overlaps.items():
                    print(ovl_start, ovl_end, count, len(features))
                    count /= len(features)
                    for p in range(ovl_start, ovl_end + 1):
                        ref_coverage[region][p] += count
                        ambig_ref_coverage[region][p] += count
                print("refcov", region, ref_coverage[region])

        for rid, regions in self.ambig_coverage.items():
            ref, reflen = bam.get_reference(rid[0] if isinstance(rid, tuple) else rid)
            for (start, end, _), overlaps in regions.items():
                features = self.db.db.get(ref, dict()).get((start, end), list())
                region = (rid, start, end)
                region_coverage = Counter({p: 0 for p in range(start, end + 1)})
                ambig_ref_coverage.setdefault(region, Counter()).update(region_coverage)
                pos_feature_dict.setdefault(region, set()).update(features)
                print("REGION", region)
                print(overlaps)
                for (ovl_start, ovl_end), count in overlaps.items():
                    for p in range(ovl_start, ovl_end + 1):
                        ambig_ref_coverage[(rid, start, end)][p] += count / len(
                            features
                        )

        print("REF_COV", ref_coverage)
        print("AMBIG_REF_COV", ambig_ref_coverage)
        print("POS_FEAT", pos_feature_dict)

        domain_cov = dict()
        for pos in set(ref_coverage).union(ambig_ref_coverage):
            # features_at_position = pos_feature_dict.get(pos, set())
            for domtype in pos_feature_dict.get(pos, set()):
                domain_cov.setdefault(domtype, dict()).update(
                    {
                        "depth_uniq": list(),
                        "depth_ambig": list(),
                        "cov_uniq": list(),
                        "cov_ambig": list(),
                    }
                )

                uniq_cov = ref_coverage.get(pos, Counter())
                ambig_cov = ambig_ref_coverage.get(pos, Counter())

                if uniq_cov:
                    domain_cov[domtype]["depth_uniq"].append(
                        sum(uniq_cov.values()) / len(uniq_cov.values())
                    )
                    domain_cov[domtype]["cov_uniq"].append(
                        sum(1 for v in uniq_cov.values() if v) / len(uniq_cov.values())
                    )
                if ambig_cov:
                    domain_cov[domtype]["depth_ambig"].append(
                        sum(ambig_cov.values()) / len(ambig_cov.values())
                    )
                    domain_cov[domtype]["cov_ambig"].append(
                        sum(1 for v in ambig_cov.values() if v) / len(ambig_cov.values())
                    )

        with gzip.open(self.out_prefix + ".covsum.txt.gz", "wt") as cov_out:
            print(
                "#domain",
                "depth_unique",
                "depth_combined",
                "coverage_unique",
                "coverage_combined",
                sep="\t",
                file=cov_out,
            )
            for domtype, counts in sorted(domain_cov.items()):

                depth_uniq, depth_ambig = [
                    c for c in counts.get("depth_uniq", list()) if c is not None
                ], [c for c in counts.get("depth_ambig", list()) if c is not None]

                depth_uniq_ = (sum(depth_uniq) / len(depth_uniq)) if depth_uniq else 0
                print(
                    domtype,
                    "UNIQ",
                    depth_uniq,
                    sum(depth_uniq),
                    len(depth_uniq),
                    "=",
                    depth_uniq_,
                )

                depth_ambig_ = (
                    (sum(depth_ambig) / len(depth_ambig)) if depth_ambig else 0
                )
                print(
                    domtype,
                    "AMBIG",
                    depth_ambig,
                    sum(depth_ambig),
                    len(depth_ambig),
                    "=",
                    depth_ambig_,
                )

                if depth_ambig_ < depth_uniq_:
                    raise ValueError(
                        f"{domtype}: depth_uniq_ + depth_ambig_ {depth_ambig_:.5f} is smaller than uniq depth {depth_uniq_:.5f}."
                    )

                cov_uniq, cov_ambig = [
                    c for c in counts.get("cov_uniq", list()) if c is not None
                ], [c for c in counts.get("cov_ambig", list()) if c is not None]

                cov_uniq_ = (sum(cov_uniq) / len(cov_uniq)) if cov_uniq else 0
                print(
                    domtype,
                    "UNIQ",
                    cov_uniq,
                    sum(cov_uniq),
                    len(cov_uniq),
                    "=",
                    cov_uniq_,
                )

                cov_ambig_ = (sum(cov_ambig) / len(cov_ambig)) if cov_ambig else 0
                print(
                    domtype,
                    "AMBIG",
                    cov_ambig,
                    sum(cov_ambig),
                    len(cov_ambig),
                    "=",
                    cov_ambig_,
                )

                if cov_ambig_ < cov_uniq_:
                    raise ValueError(
                        f"{domtype}: cov_uniq_ + cov_ambig_ {cov_ambig_:.5f} is smaller than uniq cov {cov_uniq_:.5f}."
                    )

                # print(domtype, f"{depth_uniq_:.5f}", f"{depth_ambig_:.5f}", sep="\t", file=cov_out)
                print(
                    domtype,
                    *(
                        f"{val:.5f}"
                        for val in (depth_uniq_, depth_ambig_, cov_uniq_, cov_ambig_)
                    ),
                    sep="\t",
                    file=cov_out,
                )

    def dump_counts(self, bam=None):
        print("Dumping overlap counters...", flush=True)
        print("Has ambiguous counts:", self.has_ambig_counts, flush=True)

        self._dump_feature_counts()

        if self.gene_counts:
            self._dump_gene_counts(bam=bam)

        if bam:
            self._dump_seq_counts(bam)

        if self.db.reference_type == "domain":
            self.summarise_coverage(bam)

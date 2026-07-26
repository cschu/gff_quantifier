""" module docstring """

import logging

import numpy as np


logger = logging.getLogger(__name__)


class CountMatrix:
    NUMPY_DTYPE = 'float64'  # float16 causes some overflow issue during testing
    NO_ANNOTATION = "0"  # group-tag for genes without emappper functional annotation
    RAW_COLUMNS = 0, 4
    LNORM_COLUMNS = 1, 5
    SCALED_COLUMNS = 2, 6
    RPKM_COLUMNS = 3, 7

    @classmethod
    def from_count_matrix(cls, cmatrix, rows=None):
        if rows is None:
            counts = np.array(cmatrix.counts)
            index = dict(cmatrix.index.items())
        else:
            counts = cmatrix.counts[rows, :]
            index = {}
            for (key, _), keep in zip(cmatrix.index.items(), rows):
                if keep:
                    index[key] = len(index)
            # index = {
            #     key: value
            #     for (key, value), keep in zip(cmatrix.index.items(), rows)
            #     if keep
            # }
        return cls(index=index, counts=counts)

    @staticmethod
    def calculate_scaling_factor(raw, norm):
        if norm == 0.0:
            return 1.0
        return raw / norm

    def __init__(self, ncols=2, nrows=1000, index=None, counts=None,):
        if index is not None and counts is not None:
            self.index = dict(index.items())
            self.counts = counts
        else:
            self.index = {}
            self.counts = np.zeros(
                (nrows, ncols,),
                dtype=CountMatrix.NUMPY_DTYPE,
            )

    def has_record(self, key):
        return self.index.get(key) is not None

    def _resize(self):
        nrows, ncols = self.counts.shape[0], self.counts.shape[1] 
        if len(self.index) == nrows:
            self.counts = np.pad(
                self.counts,
                ((0, nrows + 1000), (0, 0),),  # pad_width, s. numpy.pad: (0 new rows at the start, 1000 new rows at the end; no new cols on either side)
            )
        return len(self.index)

    def __getitem__(self, key):
        key_index = self.index.get(key)
        if key_index is None:
            key_index = self.index[key] = self._resize()
        return self.counts[key_index]

    def __setitem__(self, key, value):
        key_index = self.index.get(key)
        if key_index is None:
            key_index = self.index[key] = self._resize()
        self.counts[key_index] = value

    def __iter__(self):
        yield from zip(self.index.keys(), self.counts)

    def sum(self):
        return self.counts.sum(axis=0)

    def scale_column(self, col, target_col, factor, rows=None):
        # apply scaling factors
        if rows is None:
            self.counts[:, target_col] = self.counts[:, col] * factor
        else:
            self.counts[rows, target_col] = self.counts[rows, col] * factor

    def drop_unindexed(self):
        self.counts = self.counts[0:len(self.index), :]

    def dump(self, prefix="prefix", state="genes", labels=None,):
        with open(f"{prefix}.cmatrix.{state}.txt", "wt", encoding="UTF-8",) as _out:
            if labels is None:
                for index, counts in self:
                    print(index, *counts, sep="\t", file=_out)
            else:
                for (index, counts), label in zip(self, labels):
                    print(label, *counts, sep="\t", file=_out)

    def colsum(self, col):
        return self.counts[:, col].sum()

    def colsums(self):
        return self.counts.sum(axis=0)

    # def generate_gene_counts(self, lengths):
    #     logger.info("LENGTHS ARRAY = %s", lengths.shape)
    #     logger.info("INDEX SIZE = %s", len(self.index))

    #     # remove the un-indexed rows
    #     counts = self.counts[0:len(self.index), :]

    #     # calculate combined_raw
    #     counts[:, 1:2] += counts[:, 0:1]
        
    #     # duplicate the raw counts
    #     counts = np.column_stack(
    #         (
    #             counts[:, 0], counts[:, 0], counts[:, 0], counts[:, 0], # 0, 1, 2, 3
    #             counts[:, 1], counts[:, 1], counts[:, 1], counts[:, 1], # 4, 5, 6, 7
    #         ),
    #     )

    #     # length-normalise the lnorm columns
    #     counts[:, CountMatrix.LNORM_COLUMNS[0]::4] /= lengths[:, None]

    #     count_sums = counts.sum(axis=0)

    #     uniq_scaling_factor, combined_scaling_factor = (
    #         CountMatrix.calculate_scaling_factor(*count_sums[CountMatrix.RAW_COLUMNS[0]:CountMatrix.LNORM_COLUMNS[0]]),
    #         CountMatrix.calculate_scaling_factor(*count_sums[CountMatrix.RAW_COLUMNS[1]:CountMatrix.LNORM_COLUMNS[1]]),
    #     )

    #     logger.info(
    #         "AC:: TOTAL GENE COUNTS: uraw=%s unorm=%s craw=%s cnorm=%s => SF: %s %s",
    #         count_sums[0], count_sums[1], count_sums[4], count_sums[5],
    #         uniq_scaling_factor, combined_scaling_factor,
    #     )

    #     # apply scaling factors
    #     counts[:, CountMatrix.SCALED_COLUMNS[0]] = counts[:, CountMatrix.LNORM_COLUMNS[0]] * uniq_scaling_factor
    #     counts[:, CountMatrix.SCALED_COLUMNS[1]] = counts[:, CountMatrix.LNORM_COLUMNS[1]] * combined_scaling_factor

    #     self.counts = counts

    #     return self

    def to_full_count_matrix(self):
        # remove the un-indexed rows
        counts = self.counts[0:len(self.index), :]

        # calculate combined_raw
        counts[:, 1:2] += counts[:, 0:1]
                
        # duplicate the raw counts
        counts = np.column_stack(
            (
                counts[:, 0], counts[:, 0], counts[:, 0], counts[:, 0], # 0, 1, 2, 3
                counts[:, 1], counts[:, 1], counts[:, 1], counts[:, 1], # 4, 5, 6, 7
            ),
        )

        return counts


class GeneCountMatrix(CountMatrix):
    def __init__(self, m: CountMatrix, lengths=None,):
        CountMatrix.__init__(self, index=m.index, counts=m.to_full_count_matrix(),)

        if lengths is not None:
            # length-normalise the lnorm columns
            self.counts[:, CountMatrix.LNORM_COLUMNS[0]::4] /= lengths[:, None]
            
            count_sums = self.counts.sum(axis=0)
            
            uniq_scaling_factor, combined_scaling_factor = (
                CountMatrix.calculate_scaling_factor(*count_sums[CountMatrix.RAW_COLUMNS[0]:CountMatrix.LNORM_COLUMNS[0] + 1]),
                CountMatrix.calculate_scaling_factor(*count_sums[CountMatrix.RAW_COLUMNS[1]:CountMatrix.LNORM_COLUMNS[1] + 1]),
            )

            logger.info(
                "AC:: TOTAL GENE COUNTS: uraw=%s unorm=%s craw=%s cnorm=%s => SF: %s %s",
                count_sums[0], count_sums[1], count_sums[4], count_sums[5],
                uniq_scaling_factor, combined_scaling_factor,
            )

            # apply scaling factors
            self.counts[:, CountMatrix.SCALED_COLUMNS[0]] = self.counts[:, CountMatrix.LNORM_COLUMNS[0]] * uniq_scaling_factor
            self.counts[:, CountMatrix.SCALED_COLUMNS[1]] = self.counts[:, CountMatrix.LNORM_COLUMNS[1]] * combined_scaling_factor
        

    def group_gene_counts(self, ggroups):
        ggroup_counts = CountMatrix(ncols=8)
        for (_, gene_counts), ggroup_id in zip(self, ggroups):
            ggroup_counts[ggroup_id] += gene_counts

        return ggroup_counts

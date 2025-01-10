""" module docstring """

import logging

import numpy as np


logger = logging.getLogger(__name__)


class CountMatrix:
    NUMPY_DTYPE = 'float64'  # float16 causes some overflow issue during testing

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
        nrows = self.counts.shape[0]
        if len(self.index) == nrows:
            self.counts = np.pad(
                self.counts,
                ((0, nrows + 1000), (0, 0),),
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

    def scale_column(self, col_index, factor, rows=None):
        # apply scaling factors
        if rows is None:
            self.counts[:, col_index + 1] = self.counts[:, col_index] * factor
        else:
            self.counts[rows, col_index + 1] = self.counts[rows, col_index] * factor

    def drop_unindexed(self):
        self.counts = self.counts[0:len(self.index), :]

    def generate_gene_counts(self, lengths):
        logger.info("LENGTHS ARRAY = %s", lengths.shape)
        logger.info("INDEX SIZE = %s", len(self.index))

        # remove the un-indexed rows
        counts = self.counts[0:len(self.index), :]

        # calculate combined_raw
        counts[:, 1:2] += counts[:, 0:1]

        # duplicate the raw counts
        counts = np.column_stack(
            (
                counts[:, 0], counts[:, 0], counts[:, 0],  # 0, 1, 2
                counts[:, 1], counts[:, 1], counts[:, 1],  # 3, 4, 5
            ),
        )

        # length-normalise the lnorm columns
        counts[:, 1::3] /= lengths[:, None]

        count_sums = counts.sum(axis=0)

        uniq_scaling_factor, combined_scaling_factor = (
            CountMatrix.calculate_scaling_factor(*count_sums[0:2]),
            CountMatrix.calculate_scaling_factor(*count_sums[3:5]),
        )

        logger.info(
            "AC:: TOTAL GENE COUNTS: uraw=%s unorm=%s craw=%s cnorm=%s => SF: %s %s",
            count_sums[0], count_sums[1], count_sums[3], count_sums[4],
            uniq_scaling_factor, combined_scaling_factor,
        )

        # apply scaling factors
        counts[:, 2] = counts[:, 1] * uniq_scaling_factor
        counts[:, 5] = counts[:, 4] * combined_scaling_factor

        self.counts = counts

        return self
    
    def dump(self, state="genes", labels=None,):
        with open(f"CountMatrix.{state}.txt", "wt") as _out:
            if labels is None:
                for index, counts in self:
                    print(index, *counts, sep="\t", file=_out)
            else:
                for (index, counts), label in zip(self, labels):
                    print(label, *counts, sep="\t", file=_out)


    def group_gene_counts(self, ggroups):
        ggroup_index = {}
        # for gene_id, gene_counts in self:
        #     ggroup_id = gene_id.split(".")[-1]
        #     g_key_index = ggroup_index.get(ggroup_id)
        for (_, gene_counts), ggroup_id in zip(self, ggroups):
            g_key_index = ggroup_index.get(ggroup_id)
            # gene_counts = self.counts[self.index[key]]
            if g_key_index is None:
                g_key_index = ggroup_index[ggroup_id] = len(ggroup_index)
                self.counts[g_key_index] = gene_counts
                # logger.info("CM.group_gene_counts: Adding %s to new group %s (%s).", str(gene_counts), ggroup_id, g_key_index)
            else:
                self.counts[g_key_index] += gene_counts
                # logger.info("CM.group_gene_counts: Adding %s to group %s (%s).", str(gene_counts), ggroup_id, g_key_index)

        # replace index with grouped index
        self.index = ggroup_index

        # remove the un-indexed (ungrouped) rows
        self.counts = self.counts[0:len(self.index), :]

        return self

    def colsum(self, col):
        return self.counts[:, col].sum()

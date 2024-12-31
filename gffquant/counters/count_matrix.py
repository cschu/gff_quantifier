import logging

import numpy as np


logger = logging.getLogger(__name__)

class CountMatrix:

    @staticmethod
    def calculate_scaling_factor(raw, norm):
        if norm == 0.0:
            return 1.0
        return raw / norm

    def __init__(self, ncols, nrows=1000):
        self.index = {}
        self.counts = np.zeros(
            (nrows, ncols,),
            dtype='float64',
        )

    def _resize(self):
        nrows = self.counts.shape[0]
        if len(self.index) == nrows:
            self.counts = np.pad(
                self.counts,
                ((0, nrows * 2), (0, 0),),
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
    
    def group_gene_counts(self, ggroups):
        ggroup_index = {}
        for (key, key_index), ggroup_id in zip(self.index.items(), ggroups):
            g_key_index = ggroup_index.get(ggroup_id)
            gene_counts = self.counts[self.index[key]]
            if g_key_index is None:
                g_key_index = ggroup_index[ggroup_id] = len(ggroup_index)
                self.counts[g_key_index] = gene_counts
            else:
                self.counts[g_key_index] += gene_counts

        # replace index with grouped index
        self.index = ggroup_index

        # remove the un-indexed (ungrouped) rows
        self.counts = self.counts[0:len(self.index), :]

        return self
    
    def colsum(self, col):
        return self.counts[:, col].sum()
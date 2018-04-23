import os

from definitions import ROOT_DIR

import numpy as np
import pandas as pd
import networkx as nx
import dask.dataframe as dd
from dask.multiprocessing import get

from genomic import GenomicData


class Association():
    def __init__(self, A: GenomicData, B: GenomicData, C, bi_direction=False):
        self.A = A
        self.B = B
        self.C = C

        self.W = nx.DiGraph()

    def fit(self, putative_assocs, map_function, n_jobs=4):
        edges_added = 0

        if putative_assocs is not None:
            putative_dd = dd.from_pandas(putative_assocs, npartitions=n_jobs)

            res = putative_dd.map_partitions(map_function, meta=putative_dd).compute(get=get)

            for res_partition in res: 
                for tup in res_partition:
                    self.W.add_edge(tup[0], tup[1], dys=tup[2])
                    edges_added += 1

        return edges_added

    def calc_dys_A_B(df, miRNA_A, miRNA_B, gene_A, gene_B):
        result = []
        for row in df.iterrows():
            m = row[1]['MiRBase ID']
            t = row[1]['Gene Symbol']
            miRNA_gene_A_corr = np.dot(miRNA_A[m] - np.mean(miRNA_A[m]),
                                       gene_A[t] - np.mean(gene_A[t])) / \
                                ((n_A - 1) * np.std(miRNA_A[m]) * np.std(gene_A[t]))
            miRNA_gene_B_corr = np.dot(miRNA_B[m] - np.mean(miRNA_B[m]),
                                       gene_B[t] - np.mean(gene_B[t])) / \
                                ((n_B - 1) * np.std(miRNA_B[m]) * np.std(gene_B[t]))
            dys = miRNA_gene_A_corr - miRNA_gene_B_corr
            p_value = self.z_to_p_value(self.fisher_r_to_z(miRNA_gene_A_corr, n_A, miRNA_gene_B_corr, n_B))

            if p_value <= p_threshold and (miRNA_gene_A_corr < 0 or miRNA_gene_B_corr < 0):
                result.append((m, t, p_value))
        return result

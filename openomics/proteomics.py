import networkx as nx
import pandas as pd

from .database.annotation import Annotatable
from .transcriptomics import ExpressionData


class ProteinExpression(ExpressionData, Annotatable):
    def __init__(self, cohort_name, file_path, columns, genes_col_name, gene_index, sample_index="sample_barcode",
                 transposed=True,
                 log2_transform=False, npartitions=0):
        super().__init__(cohort_name, file_path, columns=columns, genes_col_name=genes_col_name, gene_index=gene_index,
                         sample_index=sample_index, transposed=transposed, log2_transform=log2_transform,
                         npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__

    def process_HPRD_PPI_network(self, ppi_data_file_path):
        HPRD_PPI = pd.read_table(ppi_data_file_path, header=None)
        self.HPRD_PPI_network = nx.from_pandas_edgelist(HPRD_PPI, source=0, target=3,
                                                        create_using=nx.DiGraph())

    def get_HPRD_PPI_network_edgelist(self):
        return self.HPRD_PPI_network.edges(data=True)

    def process_STRING_PPI_network(self, ppi_data_file_path):
        pass

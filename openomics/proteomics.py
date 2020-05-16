import networkx as nx
import pandas as pd

from .database import Annotatable
from .transcriptomics import ExpressionData


class Protein(ExpressionData, Annotatable):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None,
                 sample_index_by="sample_index", log2_transform=False, npartitions=None):
        super(Protein, self).__init__(cohort_name, data, transposed=transposed, columns=columns,
                                      gene_index_by=gene_index_by,
                                      sample_index_by=sample_index_by, log2_transform=log2_transform,
                                      npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__

    def process_HPRD_PPI_network(self, ppi_data_file_path):
        HPRD_PPI = pd.read_table(ppi_data_file_path, header=None)
        self.HPRD_PPI_network = nx.from_pandas_edgelist(HPRD_PPI, source=0, target=3,
                                                        create_using=nx.DiGraph())


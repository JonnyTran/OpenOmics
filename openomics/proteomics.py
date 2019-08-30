import networkx as nx

from openomics.database.annotation import Annotatable
from openomics.transcriptomics import ExpressionData


class Protein(ExpressionData, Annotatable):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super().__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                         transposed=transposed, log2_transform=log2_transform)

    @classmethod
    def name(self):
        return "PRO"

    def process_HPRD_PPI_network(self, ppi_data_file_path):
        HPRD_PPI = pd.read_table(ppi_data_file_path, header=None)
        self.HPRD_PPI_network = nx.from_pandas_edgelist(HPRD_PPI, source=0, target=3,
                                                        create_using=nx.DiGraph())

    def get_HPRD_PPI_network_edgelist(self):
        return self.HPRD_PPI_network.edges(data=True)

    def process_STRING_PPI_network(self, ppi_data_file_path):
        pass

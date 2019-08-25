import os
import pandas as pd
import networkx as nx
from openTCGA.annotation import *


class Interactions(Database):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_interactions(self, *args):
        pass


class LncBase(Interactions, Database):
    def __init__(self, import_folder) -> None:
        if not os.path.isdir(import_folder) or not os.path.exists(import_folder):
            raise NotADirectoryError(import_folder)
        self.folder_path = import_folder

        self.df = pd.read_table(os.path.join(import_folder, "lncBase", "lncBaseV2_predicted_human_data.csv"))


    def get_rename_dict(self, from_index, to_index):
        lncBase_gene_id_to_name_dict = pd.Series(self.df["geneName"].values,
                                                 index=self.df["geneId"]).to_dict()
        return lncBase_gene_id_to_name_dict

    def get_interactions(self, organism="Homo sapiens", tissue=None, rename_dict=None, data=True):
        lncbase_df = self.df

        lncbase_df = lncbase_df[lncbase_df["species"] == organism]
        if tissue is not None:
            lncbase_df = lncbase_df[lncbase_df["tissue"] == tissue]

        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(lncbase_df, source='mirna', target='geneId',
                                                               edge_attr=["tissue", "positive_negative"],
                                                               create_using=nx.DiGraph())

        if rename_dict is not None:
            lncBase_lncRNA_miRNA_network = nx.relabel_nodes(lncBase_lncRNA_miRNA_network, rename_dict)

        return lncBase_lncRNA_miRNA_network.edges(data=data)
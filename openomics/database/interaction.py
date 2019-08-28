import networkx as nx

from openomics.database.annotation import *


class Interactions(Database):
    __metaclass__ = ABCMeta

    def get_network_edgelist(self):
        if hasattr(self, "network"):
            return self.network.edges(data=True)
        else:
            print(self.__class__.__str__(), "does not have network interaction data yet. (at self.network)")
            return None


    @abstractmethod
    def get_interactions(self, source_index, target_index, edge_attr,  *args): raise NotImplementedError


class LncBase(Interactions, Database):
    def __init__(self, import_folder, file_resources, column_rename_dict=None) -> None:
        if file_resources is None:
            file_resources = {}
            file_resources["LncBasev2_download.csv"] = os.path.join(import_folder, "LncBasev2_download.csv")

        super().__init__(import_folder, file_resources, column_rename_dict)

    def load_data(self, file_resources) -> pd.DataFrame:
        return pd.read_table(file_resources["LncBasev2_download.csv"], low_memory=True)

    def get_rename_dict(self, from_index, to_index):
        gene_id_to_gene_name_dict = pd.Series(self.df["geneName"].values,
                                                 index=self.df["geneId"]).to_dict()
        return gene_id_to_gene_name_dict

    def get_interactions(self, source_index="mirna", target_index="geneId", edge_attr=["tissue", "positive_negative"],
                         organism="Homo sapiens", tissue=None, rename_dict=None, ):
        lncbase_df = self.df

        lncbase_df = lncbase_df[lncbase_df["species"] == organism]
        if tissue is not None:
            lncbase_df = lncbase_df[lncbase_df["tissue"] == tissue]

        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(lncbase_df, source=source_index, target=target_index,
                                                               edge_attr=edge_attr,
                                                               create_using=nx.DiGraph())
        if rename_dict is not None:
            lncBase_lncRNA_miRNA_network = nx.relabel_nodes(lncBase_lncRNA_miRNA_network, rename_dict)

        if edge_attr is None:
            data = False
        else:
            data = True

        return lncBase_lncRNA_miRNA_network.edges(data=data)


class lncRInter(Interactions, Database):
    pass

class LncRNATarget(Interactions, Database):
    pass

class lncRNome(Interactions, Database):
    pass


class NPInter(Interactions, Database):
    pass


class TargetScan(Interactions, Database):
    pass



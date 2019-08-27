import networkx as nx

from openTCGA.database.annotation import *


class Interactions(Database):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_interactions(self, source_index, target_index, edge_attr,  *args): raise NotImplementedError


class LncBase(Interactions, Database):
    def __init__(self, import_folder, file_resources, column_rename_dict) -> None:
        if not os.path.isdir(import_folder) or not os.path.exists(import_folder):
            raise NotADirectoryError(import_folder)
        self.folder_path = import_folder
        self.df = pd.read_table(os.path.join(import_folder, "LncBasev2_download.csv"))
        print(self.df.columns.tolist())

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



import networkx as nx

from openomics.database.annotation import *


class Interactions(Dataset):
    __metaclass__ = ABCMeta

    def __init__(self, import_folder, file_resources, source_index, target_index, edge_attr=None, rename_dict=None,
                 **kwargs):
        """
        This is an abstract class used to instantiate a database given a folder containing various file resources. When creating a Database class, the load_data function is called where the file resources are load as a DataFrame and performs necessary processings. This class provides an interface for RNA classes to annotate various genomic annotations, functional annotations, sequences, and disease associations.
        Args:
            import_folder (str):
                The folder path containing the data files
            file_resources (dict): default None,
                Used to list required files for load_network of the dataset. A dictionary where keys are required filenames and value are file paths. If None, then the class constructor should automatically build the required file resources dict.
            source_index (str): Column name of DataFrame to be used as the source node names
            target_index (str): Column name of DataFrame to be used as the target node names
            edge_attr (list): A list of attributes to
            col_rename (dict): default None,
                A dictionary to rename columns in the data table. If None, then automatically load defaults.
            **kwargs: Additional arguments that may be passed to load_data function
        """
        if not os.path.isdir(import_folder) or not os.path.exists(import_folder):
            raise NotADirectoryError(import_folder)
        else:
            for _, filepath in file_resources.items():
                if not os.path.exists(filepath):
                    raise FileNotFoundError(filepath)

        self.import_folder = import_folder
        self.file_resources = file_resources
        self.network = self.load_network(file_resources, source_index, target_index, edge_attr, **kwargs)
        if self.network is None:
            raise Exception("Make sure load_network() returns a Networkx Graph")
        if rename_dict is not None:
            self.network = nx.relabel_nodes(self.network, rename_dict)
        print("{}: {}".format(self.name(), nx.info(self.network)))

    @abstractmethod
    def load_network(self, file_resources, source_index, target_index, edge_attr, **kwargs) -> nx.Graph:
        raise NotImplementedError

    def get_interactions(self, edge_attr=False):
        if hasattr(self, "network"):
            return self.network.edges(data=edge_attr)
        else:
            raise Exception(
                "{} does not have network interaction data yet. Must run load_network() first.".format(self.name()))


class LncBase(Interactions, Dataset):
    def __init__(self, import_folder, file_resources=None, source_index="mirna", target_index="geneId",
                 edge_attr=["tissue", "positive_negative"],
                 rename_dict=None, organism="Homo sapiens", tissue=None) -> None:
        """

        Args:
            import_folder (str):
            file_resources (dict): default None.
            col_rename (dict): default None.
            species (str): {'Homo sapiens', "Kaposi's sarcoma-associated herpesvirus (KSHV)", 'Epstein–Barr virus', 'Mus musculus'}
        """
        if file_resources is None:
            file_resources = {}
            file_resources["LncBasev2_download.csv"] = os.path.join(import_folder, "LncBasev2_download.csv")

        super().__init__(import_folder, file_resources, source_index, target_index, edge_attr, rename_dict,
                         organism=organism, tissue=tissue)

    def get_rename_dict(self, from_index="geneId", to_index="geneName"):
        lncbase_df = pd.read_table(self.file_resources["LncBasev2_download.csv"], low_memory=True)
        gene_id_to_gene_name_dict = pd.Series(lncbase_df["geneName"].values,
                                              index=lncbase_df["geneId"]).to_dict()
        return gene_id_to_gene_name_dict

    def load_network(self, file_resources, source_index="mirna", target_index="gene_id",
                     edge_attr=["tissue", "positive_negative"],
                     organism="Homo sapiens", tissue=None):
        lncbase_df = pd.read_table(file_resources["LncBasev2_download.csv"], low_memory=True)
        lncbase_df.replace({"species": {"Homo Sapiens": "Homo sapiens", "Mus Musculus": "Mus musculus"}}, inplace=True)

        lncbase_df = lncbase_df[lncbase_df["species"] == organism]
        if tissue is not None:
            lncbase_df = lncbase_df[lncbase_df["tissue"] == tissue]

        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(lncbase_df, source=source_index, target=target_index,
                                                               edge_attr=edge_attr,
                                                               create_using=nx.DiGraph())
        return lncBase_lncRNA_miRNA_network


class lncRInter(Interactions, Dataset):
    pass


class LncRNATarget(Interactions, Dataset):
    pass


class lncRNome(Interactions, Dataset):
    pass


class NPInter(Interactions, Dataset):
    pass


class MiRTarBase(Interactions):
    def __init__(self, import_folder, file_resources=None, source_index="miRNA", target_index="Target Gene",
                 edge_attr=["Support Type"], rename_dict=None, species="Homo sapiens"):
        if file_resources is None:
            file_resources = {}
            file_resources["miRTarBase_MTI.xlsx"] = os.path.join(import_folder, "miRTarBase_MTI.xlsx")

        super().__init__(import_folder, file_resources, source_index, target_index, edge_attr, rename_dict,
                         species=species)

    def load_network(self, source_index="miRNA", target_index="Target Gene", edge_attr=["Support Type"],
                     rename_dict=None, species="Homo sapiens"):
        table = pd.read_excel(self.file_resources["miRTarBase_MTI.xlsx"])
        if species:
            table = table[table["Species (Target Gene)"].str.lower() == species.lower()]
        # table['miRNA'] = table['miRNA'].str.lower()
        # table['miRNA'] = table['miRNA'].str.replace("-3p.*|-5p.*", "")
        mir_target_network = nx.from_pandas_edgelist(table, source=source_index, target=target_index,
                                                     edge_attr=edge_attr,
                                                     create_using=nx.DiGraph())
        return mir_target_network


class TargetScan(Interactions, Dataset):
    def __init__(self, import_folder, file_resources=None, source_index="mirna", target_index="gene_id",
                 edge_attr=["tissue", "positive_negative"], rename_dict=None, species=9606):

        if file_resources is None:
            file_resources = {}
            file_resources["miR_Family_Info.txt"] = os.path.join(import_folder, "miR_Family_Info.txt")
            file_resources["Predicted_Targets_Info.default_predictions.txt"] = os.path.join(import_folder,
                                                                                            "Predicted_Targets_Info.default_predictions.txt")

        super().__init__(import_folder, file_resources, source_index, target_index, edge_attr, rename_dict,
                         species=species)

    def load_network(self, file_resources, source_index="mirna", target_index="gene_id",
                     edge_attr=["tissue", "positive_negative"], species=9606):
        self.df = self.process_miR_family_info(file_resources, species)
        interactions_df = self.process_targetScan_mirna_target_interactions(file_resources, self.df, species)
        mir_target_network = nx.from_pandas_edgelist(interactions_df, source="MiRBase ID", target="Gene Symbol",
                                                     create_using=nx.DiGraph())
        return mir_target_network

    def process_miR_family_info(self, file_resources, species):
        miR_Family_Info_df = pd.read_table(file_resources["miR_Family_Info.txt"], delimiter='\t')
        if species:
            miR_Family_Info_df = miR_Family_Info_df[miR_Family_Info_df['Species ID'] == species]
        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        miR_Family_Info_df['MiRBase ID'] = miR_Family_Info_df['MiRBase ID'].str.lower()
        miR_Family_Info_df['MiRBase ID'] = miR_Family_Info_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        miR_Family_Info_df.drop_duplicates(inplace=True)
        miR_Family_Info_df = miR_Family_Info_df[
            ['miR family', 'MiRBase ID', 'Seed+m8', 'Mature sequence', 'Family Conservation?', 'MiRBase Accession']]
        miR_Family_Info_df['MiRBase ID'] = miR_Family_Info_df['MiRBase ID'].astype(str)
        miR_Family_Info_df.set_index("MiRBase ID", inplace=True)
        miR_Family_Info_df.drop('MiRBase ID', axis=1, inplace=True)
        return miR_Family_Info_df

    def process_targetScan_mirna_target_interactions(self, file_resources, targetScan_family_df, species):
        # Load data frame from file
        interactions_df = pd.read_table(file_resources["Predicted_Targets_Info.default_predictions.txt"],
                                        delimiter='\t', low_memory=True)

        # Select only homo sapiens miRNA-target pairs
        if species:
            interactions_df = interactions_df[interactions_df["Species ID"] == species]
            targetScan_family_df = targetScan_family_df[targetScan_family_df['Species ID'] == species]

        interactions_df = interactions_df.filter(items=["miR Family", "Gene Symbol"], axis="columns")
        targetScan_family_df = targetScan_family_df.filter(items=['miR family', 'MiRBase ID'], axis="columns")

        # map miRBase ID names to miR Family
        targetScan_family_df.rename(columns={'miR family': 'miR Family'}, inplace=True)
        interactions_df = pd.merge(interactions_df, targetScan_family_df, how='inner', on="miR Family")
        interactions_df = interactions_df[["MiRBase ID", "Gene Symbol"]]

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        interactions_df['MiRBase ID'] = interactions_df['MiRBase ID'].str.lower()
        interactions_df['MiRBase ID'] = interactions_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")
        interactions_df.drop_duplicates(inplace=True)
        return interactions_df

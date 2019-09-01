from openomics.database.annotation import *
import networkx as nx

class Interactions(Dataset):
    __metaclass__ = ABCMeta

    def __init__(self, import_folder, file_resources, source_col_name, target_col_name, source_index, target_index,
                 edge_attr=None, directed=True, rename_dict=None, **kwargs):
        """
        This is an abstract class used to instantiate a database given a folder containing various file resources. When creating a Database class, the load_data function is called where the file resources are load as a DataFrame and performs necessary processings. This class provides an interface for RNA classes to annotate various genomic annotations, functional annotations, sequences, and disease associations.
        Args:
            import_folder (str):
                The folder path containing the data files.
            file_resources (dict): default None,
                Used to list required files for load_network of the dataset. A dictionary where keys are required filenames and value are file paths. If None, then the class constructor should automatically build the required file resources dict.
            source_col_name (str):
                Column name of DataFrame to be used as the source node names.
            target_col_name (str):
                Column name of DataFrame to be used as the target node names.
            edge_attr (list):
                A list of column names to be included as attributes for each edge (source-target pairs).
            directed (bool): default True,
                Whether to create a directed or an undirected network.
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
        self.source_index = source_index
        self.target_index = target_index
        self.network = self.load_network(file_resources=file_resources, source_col_name=source_col_name,
                                         target_col_name=target_col_name,
                                         edge_attr=edge_attr, directed=directed, **kwargs)
        self.network.name = self.name()

        if self.network is None:
            raise Exception(
                "Make sure load_network() returns a Networkx Graph and is called with super().__init__() in the constructor.")

        if rename_dict is not None:
            self.network = nx.relabel_nodes(self.network, rename_dict)

        print("{}".format(nx.info(self.network)))

    @abstractmethod
    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, **kwargs):
        raise NotImplementedError

    def get_interactions(self, nodelist, data=False, inclusive=False):
        """

        Args:
            nodelist (list):
                A list of nodes to fetch edges from
            data (bool): default False
                Whether to include edge attributes
            inclusive (bool): default False
                Whether to only retrieve edges from nodes inclusive in nodelist.

        Returns:
            edges (OutEdgeView): a NetworkX edgelist
        """
        if hasattr(self, "network"):
            if inclusive:
                return self.network.subgraph(nodelist).edges(data=data)
            else:
                return self.network.edges(nodelist=nodelist, data=data)

        else:
            raise Exception(
                "{} does not have network interaction data yet. Must run load_network() first.".format(self.name()))


class LncBase(Interactions, Dataset):
    def __init__(self, import_folder, file_resources=None, source_col_name="mirna", target_col_name="geneId",
                 source_index="transcript_name", target_index="gene_id",
                 edge_attr=["tissue", "positive_negative"], directed=True,
                 rename_dict=None, organism="Homo sapiens", tissue=None):
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

        super().__init__(import_folder=import_folder, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, directed=directed, rename_dict=rename_dict,
                         organism=organism, tissue=tissue)

    def get_rename_dict(self, from_index="geneId", to_index="geneName"):
        lncbase_df = pd.read_table(self.file_resources["LncBasev2_download.csv"], low_memory=True)
        gene_id_to_gene_name_dict = pd.Series(lncbase_df["geneName"].values,
                                              index=lncbase_df["geneId"]).to_dict()
        return gene_id_to_gene_name_dict

    def load_network(self, file_resources, source_col_name="mirna", target_col_name="gene_id",
                     edge_attr=["tissue", "positive_negative"], directed=True,
                     organism="Homo sapiens", tissue=None):
        df = pd.read_table(file_resources["LncBasev2_download.csv"], low_memory=True)
        print(self.name(), df.columns.tolist())
        df.replace({"species": {"Homo Sapiens": "Homo sapiens", "Mus Musculus": "Mus musculus"}}, inplace=True)

        if organism is not None:
            df = df[df["species"].str.lower() == organism.lower()]
        if tissue is not None:
            df = df[df["tissue"].str.lower() == tissue.lower()]

        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                                               edge_attr=edge_attr,
                                                               create_using=nx.DiGraph() if directed else nx.Graph())
        return lncBase_lncRNA_miRNA_network


class lncRInter(Interactions):
    pass


class LncRNATarget(Interactions):
    pass


class lncRNome(Interactions):
    pass


class NPInter(Interactions):
    pass


class MiRTarBase(Interactions):
    COLUMNS_RENAME_DICT = {"Species (Target Gene)": "species",
                           "Support Type": "Support_Type",
                           "Target Gene": "gene_name"}

    def __init__(self, import_folder, file_resources=None, source_col_name="miRNA", target_col_name="Target Gene",
                 source_index="transcript_name", target_index="gene_name",
                 edge_attr=["Support Type"], directed=True, rename_dict=None, species="Homo sapiens",
                 strip_mirna_name=False):
        self.strip_mirna_name = strip_mirna_name

        if file_resources is None:
            file_resources = {}
            file_resources["miRTarBase_MTI.xlsx"] = os.path.join(import_folder, "miRTarBase_MTI.xlsx")

        super().__init__(import_folder=import_folder, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, directed=directed, rename_dict=rename_dict,
                         species=species)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed=True, species=None):
        df = pd.read_excel(self.file_resources["miRTarBase_MTI.xlsx"])
        print(self.name(), df.columns.tolist())

        if species:
            df = df[df["Species (Target Gene)"].str.lower() == species.lower()]

        if self.strip_mirna_name:
            df['miRNA'] = df['miRNA'].str.lower()
            df['miRNA'] = df['miRNA'].str.replace("-3p.*|-5p.*", "")

        mir_target_network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                                     edge_attr=edge_attr,
                                                     create_using=nx.DiGraph() if directed else nx.Graph())
        return mir_target_network


class TargetScan(Interactions, Dataset):
    def __init__(self, import_folder, file_resources=None, source_col_name="MiRBase ID", target_col_name="Gene Symbol",
                 source_index="transcript_name", target_index="transcript_name",
                 edge_attr=["tissue", "positive_negative"], directed=True, rename_dict=None, species=9606,
                 strip_mirna_name=False):
        self.strip_mirna_name = strip_mirna_name
        if file_resources is None:
            file_resources = {}
            file_resources["miR_Family_Info.txt"] = os.path.join(import_folder, "miR_Family_Info.txt")
            file_resources["Predicted_Targets_Info.default_predictions.txt"] = os.path.join(import_folder,
                                                                                            "Predicted_Targets_Info.default_predictions.txt")

        super().__init__(import_folder=import_folder, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, directed=directed, rename_dict=rename_dict,
                         species=species)

    def load_network(self, file_resources, source_col_name="MiRBase ID", target_col_name="Gene Symbol",
                     edge_attr=["tissue", "positive_negative"], directed=True, species=9606):
        self.df = self.process_miR_family_info_table(file_resources, species)
        interactions_df = self.process_interactions_table(file_resources, self.df, species)
        print(self.name(), interactions_df.columns.tolist())

        mir_target_network = nx.from_pandas_edgelist(interactions_df, source=source_col_name, target=target_col_name,
                                                     create_using=nx.DiGraph() if directed else nx.Graph())
        return mir_target_network

    def process_miR_family_info_table(self, file_resources, species=None):
        miR_Family_Info_df = pd.read_table(file_resources["miR_Family_Info.txt"], delimiter='\t')

        if species:
            miR_Family_Info_df = miR_Family_Info_df[miR_Family_Info_df['Species ID'] == species]

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        if self.strip_mirna_name:
            miR_Family_Info_df['MiRBase ID'] = miR_Family_Info_df['MiRBase ID'].str.lower()
            miR_Family_Info_df['MiRBase ID'] = miR_Family_Info_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")

        miR_Family_Info_df.drop_duplicates(inplace=True)
        miR_Family_Info_df = miR_Family_Info_df.filter(items=['miR family', 'MiRBase ID', 'Seed+m8', 'Mature sequence',
                                                              'Family Conservation?', 'MiRBase Accession'],
                                                       axis="columns")
        miR_Family_Info_df['MiRBase ID'] = miR_Family_Info_df['MiRBase ID'].astype(str)
        return miR_Family_Info_df

    def process_interactions_table(self, file_resources, family_to_miR_df, species):
        """
        This functions joins the interactions data table between miR Family and targets, and
        Args:
            file_resources:
            family_to_miR_df:
            species:

        Returns:

        """
        # Load data frame from file
        family_interactions_df = pd.read_table(file_resources["Predicted_Targets_Info.default_predictions.txt"],
                                               delimiter='\t', low_memory=True)

        # Select only homo sapiens miRNA-target pairs
        if species:
            family_interactions_df = family_interactions_df[family_interactions_df["Species ID"] == species]

        family_interactions_df = family_interactions_df.filter(items=["miR Family", "Gene Symbol"], axis="columns")
        family_to_miR_df = family_to_miR_df.filter(items=['miR family', 'MiRBase ID'], axis="columns")
        family_to_miR_df.rename(columns={'miR family': 'miR Family'}, inplace=True)

        # map miRBase ID names to miR Family
        # family_interactions_df = pd.merge(family_interactions_df, family_to_miR_df, how='outer', on="miR Family")

        family_to_miR_df.set_index("miR Family", inplace=True)
        family_interactions_df.set_index("miR Family", inplace=True)
        mir_interactions_df = family_interactions_df.join(family_to_miR_df, how='outer', on="miR Family").reset_index()

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        if self.strip_mirna_name:
            mir_interactions_df['MiRBase ID'] = mir_interactions_df['MiRBase ID'].str.lower()
            mir_interactions_df['MiRBase ID'] = mir_interactions_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")

        return mir_interactions_df

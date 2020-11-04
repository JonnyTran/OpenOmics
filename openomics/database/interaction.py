import networkx as nx

from openomics.database.annotation import *
from openomics.database.base import Dataset

from abc import abstractmethod


class Interactions(Dataset):
    def __init__(self, path, file_resources, source_col_name, target_col_name, source_index, target_index,
                 edge_attr=None, directed=True, rename_dict=None, npartitions=0):
        """
        This is an abstract class used to instantiate a database given a folder containing various file resources. When creating a Database class, the load_data function is called where the file resources are load as a DataFrame and performs necessary processings. This class provides an interface for RNA classes to annotate various genomic annotations, functional annotations, sequences, and disease associations.
        Args:
            path (str):
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
            npartitions:
        """
        if not os.path.isdir(path) or not os.path.exists(path):
            raise NotADirectoryError(path)
        else:
            for _, filepath in file_resources.items():
                if not os.path.exists(filepath):
                    raise FileNotFoundError(filepath)

        self.import_folder = path
        self.file_resources = file_resources
        self.source_index = source_index
        self.target_index = target_index
        self.network = self.load_network(file_resources=file_resources, source_col_name=source_col_name,
                                         target_col_name=target_col_name,
                                         edge_attr=edge_attr, directed=directed)
        self.network.name = self.name()

        if self.network is None:
            raise Exception(
                "Make sure load_network() returns a Networkx Graph and is called with super().__init__() in the constructor.")

        if rename_dict is not None:
            self.network = nx.relabel_nodes(self.network, rename_dict)

        print("{}".format(nx.info(self.network)))

    def __repr__(self):
        return f"{self.__class__.__name__}(num_nodes={self.network.number_of_nodes()}, num_edges={self.network.number_of_edges()})"

    @abstractmethod
    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed) -> nx.Graph:
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
                "{} does not have network interaction data yet. Must run load_network() and assign self.network field first.".format(
                    self.name()))


class GeneMania(Interactions):

    def __init__(self, path, file_resources=None, source_col_name="Gene_A", target_col_name="Gene_B",
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=None, directed=True, rename_dict=None):
        if edge_attr is None:
            edge_attr = ["Weight"]
        if file_resources is None:
            file_resources = {}
            file_resources["COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt"] = os.path.join(path,
                                                                                        "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt")
            file_resources["identifier_mappings.txt"] = os.path.join(path,
                                                                     "identifier_mappings.txt")

        super().__init__(path, file_resources, source_col_name, target_col_name, source_index, target_index,
                         edge_attr, directed, rename_dict)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed):
        interactions = pd.read_table(file_resources["COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt"], low_memory=True)
        identifier = pd.read_table(file_resources["identifier_mappings.txt"])

        # Rename ENSG ID's to gene names
        identifier = identifier[identifier["Source"] == "Gene Name"]
        identifier_map = pd.Series(identifier["Name"].values, index=identifier["Preferred_Name"]).to_dict()
        interactions.replace(identifier_map, inplace=True)

        genemania_RNA_RNA_network = nx.from_pandas_edgelist(interactions, source=source_col_name,
                                                            target=target_col_name,
                                                            edge_attr=edge_attr,
                                                            create_using=nx.DiGraph())
        return genemania_RNA_RNA_network


class BioGRID(Interactions):

    def __init__(self, path, file_resources=None, source_col_name="Official Symbol Interactor A",
                 target_col_name="Official Symbol Interactor B",
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=None,
                 directed=False, rename_dict=None):
        if edge_attr is None:
            edge_attr = ['Score', 'Throughput', 'Qualifications', 'Modification', 'Phenotypes']
        if file_resources is None:
            file_resources = {}
            file_resources["BIOGRID-ALL-X.X.XXX.tab2.txt"] = os.path.join(path, "BIOGRID-ALL-3.4.162.tab2.txt")

        super().__init__(path, file_resources, source_col_name, target_col_name, source_index, target_index,
                         edge_attr, directed, rename_dict)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, species=9606):
        biogrid_df = pd.read_table(file_resources["BIOGRID-ALL-X.X.XXX.tab2.txt"],
                                   na_values=["-"],
                                   usecols=['Official Symbol Interactor A',
                                            'Official Symbol Interactor B', 'Organism Interactor A', 'Score',
                                            'Throughput', 'Qualifications', 'Modification', 'Phenotypes'],
                                   low_memory=True)

        biogrid_df = biogrid_df[biogrid_df["Organism Interactor A"] == species]
        # biogrid_df = biogrid_df[biogrid_df["Throughput"] == "High Throughput"]

        biogrid_grn = nx.from_pandas_edgelist(biogrid_df, source=source_col_name, target=target_col_name,
                                              edge_attr=edge_attr,
                                              create_using=nx.DiGraph() if directed else nx.Graph())
        return biogrid_grn

class LncBase(Interactions, Dataset):
    def __init__(self, path, file_resources=None, source_col_name="mirna", target_col_name="geneId",
                 source_index="transcript_name", target_index="gene_id",
                 edge_attr=None, directed=True,
                 rename_dict=None, organism="Homo sapiens", tissue=None):
        """

        Args:
            path (str):
            file_resources (dict): default None.
        """
        self.organism = organism
        self.tissue = tissue

        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        if file_resources is None:
            file_resources = {}
            file_resources["LncBasev2_download.csv"] = os.path.join(path, "LncBasev2_download.csv")

        super(LncBase, self).__init__(path=path, file_resources=file_resources,
                                      source_col_name=source_col_name,
                                      target_col_name=target_col_name, source_index=source_index,
                                      target_index=target_index,
                                      edge_attr=edge_attr, directed=directed, rename_dict=rename_dict)

    def get_rename_dict(self, from_index="geneId", to_index="geneName"):
        lncbase_df = pd.read_table(self.file_resources["LncBasev2_download.csv"], low_memory=True)
        gene_id_to_gene_name_dict = pd.Series(lncbase_df["geneName"].values,
                                              index=lncbase_df["geneId"]).to_dict()
        return gene_id_to_gene_name_dict

    def load_network(self, file_resources, source_col_name="mirna", target_col_name="gene_id",
                     edge_attr=None, directed=True, ):
        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        df = pd.read_table(file_resources["LncBasev2_download.csv"], low_memory=True)
        print(self.name(), df.columns.tolist())
        df.replace({"species": {"Homo Sapiens": "Homo sapiens", "Mus Musculus": "Mus musculus"}}, inplace=True)

        if self.organism is not None:
            df = df[df["species"].str.lower() == self.organism.lower()]
        if self.tissue is not None:
            df = df[df["tissue"].str.lower() == self.tissue.lower()]

        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                                               edge_attr=edge_attr,
                                                               create_using=nx.DiGraph() if directed else nx.Graph())
        return lncBase_lncRNA_miRNA_network


class lncRInter(Interactions):

    def __init__(self, path, file_resources=None, source_col_name="lncrna",
                 target_col_name='Interacting partner',
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=None,
                 directed=True, rename_dict=None, organism="Homo sapiens"):
        self.organism = organism

        if edge_attr is None:
            edge_attr = ["Interaction Class", "Interaction Mode", "Tissue", "Phenotype"]
        if file_resources is None:
            file_resources = {}
            file_resources["human_interactions.txt"] = os.path.join(path, "human_interactions.txt")

        super().__init__(path, file_resources, source_col_name, target_col_name, source_index, target_index,
                         edge_attr, directed, rename_dict, )

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed):
        lncRInter_df = pd.read_table(file_resources["human_interactions.txt"])
        lncRInter_df = lncRInter_df[lncRInter_df["Organism"] == self.organism]

        # Data cleaning
        lncRInter_df.loc[lncRInter_df["Interacting partner"].str.contains("MIR"), "Interacting partner"] = \
            lncRInter_df.loc[
                lncRInter_df["Interacting partner"].str.contains("MIR"), "Interacting partner"].str.lower()
        lncRInter_df["Interacting partner"] = lncRInter_df["Interacting partner"].str.replace("mirlet", "hsa-let-")
        lncRInter_df["Interacting partner"] = lncRInter_df["Interacting partner"].str.replace("mir", "hsa-mir-")
        lncRInter_df["Interacting partner"][
            lncRInter_df["Interacting partner"].str.contains(r"[mir|let]\-[\d]+[a-z]+[\d]+")] = \
            lncRInter_df["Interacting partner"][
                lncRInter_df["Interacting partner"].str.contains(r"[mir|let]\-[\d]+[a-z]+[\d]+")].apply(
                lambda x: x[:-1] + "-" + x[-1])

        lncRInter_network = nx.from_pandas_edgelist(lncRInter_df, source=source_col_name,
                                                    target=target_col_name,
                                                    edge_attr=edge_attr,
                                                    create_using=nx.DiGraph() if directed else nx.Graph())
        return lncRInter_network


class LncRNA2Target(Interactions):
    def __init__(self, path, file_resources=None, source_col_name="lncrna_symbol",
                 target_col_name="gene_symbol",
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=None, directed=True, rename_dict=None, version="high_throughput", species=9606):
        """

        Args:
            version (str): one of ["high_throughput", "low_throughput"].
                The high_throughput version of lncRNA2Target database is v2.0 and low_throughput is v1.0, according to the database's website.
            species (str, int): one of [9606, "Homo sapiens"].
                The species column in high_throughput is formatted in int (e.g. 9606) and in low_throughput is in str (e.g. "Homo sapiens")
        """
        if edge_attr is None:
            edge_attr = ["P_Value", "direction"]
        self.version = version
        self.species = species
        if file_resources is None:
            file_resources = {}
            file_resources["lncRNA_target_from_high_throughput_experiments.txt"] = os.path.join(path,
                                                                                                "lncRNA_target_from_high_throughput_experiments.txt")

        super().__init__(path, file_resources, source_col_name, target_col_name, source_index, target_index,
                         edge_attr, directed, rename_dict, species)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed):
        if self.version == "high_throughput":
            return self.load_network_high_throughput(self, file_resources, source_col_name, target_col_name, edge_attr,
                                                     directed)
        elif self.version == "low_throughput":
            return self.load_network_low_throughput(self, file_resources, source_col_name, target_col_name, edge_attr,
                                                    directed)
        else:
            raise Exception("LncRNA2Target version argument must be one of 'high_throughput' or 'low_throughput'")

    def load_network_high_throughput(self, file_resources, source_col_name, target_col_name, edge_attr, directed):
        table = pd.read_table(file_resources["lncRNA_target_from_high_throughput_experiments.txt"], low_memory=True)
        table = table[table["species_id"] == self.species]

        table["lncrna_symbol"] = table["lncrna_symbol"].str.upper().replace("LINC", "")
        table["gene_symbol"] = table["gene_symbol"].str.upper()
        lncrna2target_high_throughput_network = nx.from_pandas_edgelist(table,
                                                                        source=source_col_name,
                                                                        target=target_col_name,
                                                                        edge_attr=edge_attr,
                                                                        create_using=nx.DiGraph() if directed else nx.Graph())
        return lncrna2target_high_throughput_network

    def load_network_low_throughput(self, file_resources, source_col_name="GENCODE_gene_name",
                                    target_col_name="Target_official_symbol",
                                    edge_attr=None, directed=True):
        table = pd.read_excel(file_resources["lncRNA_target_from_low_throughput_experiments.xlsx"])
        table = table[table["Species"] == self.species]

        table["Target_official_symbol"] = table["Target_official_symbol"].str.replace("(?i)(mir)", "hsa-mir-")
        table["Target_official_symbol"] = table["Target_official_symbol"].str.replace("--", "-")
        table["Target_official_symbol"].apply(lambda x: x.lower() if "mir" in x.lower() else x.upper())
        table["GENCODE_gene_name"] = table["GENCODE_gene_name"].str.upper()
        lncrna2target_low_throughput_network = nx.from_pandas_edgelist(table,
                                                                       source=source_col_name,
                                                                       target=target_col_name,
                                                                       edge_attr=edge_attr,
                                                                       create_using=nx.DiGraph() if directed else nx.Graph())
        return lncrna2target_low_throughput_network


class lncRNome(Interactions):
    pass


class NPInter(Interactions):
    pass


class MiRTarBase(Interactions):
    COLUMNS_RENAME_DICT = {"Species (Target Gene)": "species",
                           "Support Type": "Support_Type",
                           "Target Gene": "gene_name"}

    def __init__(self, path, file_resources=None, source_col_name="miRNA", target_col_name="Target Gene",
                 source_index="transcript_name", target_index="gene_name",
                 edge_attr=None, directed=True, rename_dict=None, species="Homo sapiens",
                 strip_mirna_name=False):
        if edge_attr is None:
            edge_attr = ["Support Type"]
        self.strip_mirna_name = strip_mirna_name
        self.species = species

        if file_resources is None:
            file_resources = {}
            file_resources["miRTarBase_MTI.xlsx"] = os.path.join(path, "miRTarBase_MTI.xlsx")

        super(MiRTarBase, self).__init__(path=path, file_resources=file_resources,
                                         source_col_name=source_col_name,
                                         target_col_name=target_col_name, source_index=source_index,
                                         target_index=target_index,
                                         edge_attr=edge_attr, directed=directed, rename_dict=rename_dict, )

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed=True):
        df = pd.read_excel(self.file_resources["miRTarBase_MTI.xlsx"])
        print(df.info())
        if self.species:
            df = df[df["Species (Target Gene)"].str.lower() == self.species.lower()]

        if self.strip_mirna_name:
            df['miRNA'] = df['miRNA'].str.lower()
            df['miRNA'] = df['miRNA'].str.replace("-3p.*|-5p.*", "")

        mir_target_network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                                     edge_attr=edge_attr,
                                                     create_using=nx.DiGraph() if directed else nx.Graph())
        return mir_target_network


class TargetScan(Interactions, Dataset):
    def __init__(self, path, file_resources=None, source_col_name="MiRBase ID", target_col_name="Gene Symbol",
                 source_index="transcript_name", target_index="transcript_name",
                 edge_attr=None, directed=True, rename_dict=None, species=9606,
                 strip_mirna_name=False):
        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        self.strip_mirna_name = strip_mirna_name
        self.species = species
        if file_resources is None:
            file_resources = {}
            file_resources["miR_Family_Info.txt"] = os.path.join(path, "miR_Family_Info.txt")
            file_resources["Predicted_Targets_Info.default_predictions.txt"] = os.path.join(path,
                                                                                            "Predicted_Targets_Info.default_predictions.txt")

        super(TargetScan, self).__init__(path=path, file_resources=file_resources,
                                         source_col_name=source_col_name,
                                         target_col_name=target_col_name, source_index=source_index,
                                         target_index=target_index,
                                         edge_attr=edge_attr, directed=directed, rename_dict=rename_dict)

    def load_network(self, file_resources, source_col_name="MiRBase ID", target_col_name="Gene Symbol",
                     edge_attr=None, directed=True):
        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        self.df = self.process_miR_family_info_table(file_resources, self.species)
        interactions_df = self.process_interactions_table(file_resources, self.df, self.species)
        print(self.name(), interactions_df.columns.tolist())

        mir_target_network = nx.from_pandas_edgelist(interactions_df,
                                                     source=source_col_name, target=target_col_name,
                                                     edge_attr=edge_attr,
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

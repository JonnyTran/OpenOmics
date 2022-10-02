import copy
import warnings
from abc import abstractmethod
from collections.abc import Iterable
from typing import List, Dict, Any, Union

import networkx as nx
from Bio import SeqIO
from logzero import logger

from openomics.database.annotation import *
from openomics.database.base import Database
from openomics.database.sequence import SequenceDatabase
from openomics.transforms.df import filter_rows


class Interactions(Database):
    def __init__(self, path, file_resources: Dict, source_col_name: str = None, target_col_name: str = None,
                 source_index: str = None, target_index: str = None, edge_attr: List[str] = None, filters: dict = None,
                 directed: bool = True, relabel_nodes: dict = None, blocksize=None, **kwargs):
        """
        This is an abstract class used to instantiate a database given a folder containing various file resources. When creating a Database class, the load_data function is called where the file resources are load as a DataFrame and performs necessary processings. This class provides an interface for RNA classes to annotate various genomic annotation, functional annotation, sequences, and disease associations.
        Args:
            path (str):
                The folder path containing the data files.
            file_resources (dict):
                Default None, used to list required files for load_network of the dataset. A dictionary where keys are required filenames and value are file paths. If None, then the class constructor should automatically build the required file resources dict.
            source_col_name (str):
                Column name of DataFrame to be used as the source node names.
            target_col_name (str):
                Column name of DataFrame to be used as the target node names.
            source_index (str):
                One of {"gene_name", "gene_id", "transcript_name", "transcript_id", "protein_name", "protein_id"}
            target_index (str):
                One of {"gene_name", "gene_id", "transcript_name", "transcript_id", "protein_name", "protein_id"}
            edge_attr (list):
                A list of column names to be included as attributes for each edge (source-target pairs).
            filters (dict):
                Optional. A dict with key matching the data table (from load_network()) columns and values for the filtering on that column.
            directed (bool): default True,
                Whether to create a directed or an undirected network.
            relabel_nodes (dict): default None,
                A dictionary to rename nodes in the network, where the nodes with name <dict[key]> will be renamed to <dict[value]>
            blocksize ():
        """
        # This class should NOT call super's __init__()
        super().__init__(path=path, file_resources=file_resources, blocksize=blocksize, **kwargs)
        self.source_index = source_index
        self.target_index = target_index
        self.network = self.load_network(file_resources=self.file_resources, source_col_name=source_col_name,
                                         target_col_name=target_col_name, edge_attr=edge_attr, directed=directed,
                                         filters=filters, blocksize=blocksize)

        # self.network.name = self.name()
        # if self.network is None:
        #     raise Exception(
        #         "Make sure load_network() returns a Networkx Graph and is called with super().__init__() in the constructor.")

        if relabel_nodes is not None:
            self.network = nx.relabel_nodes(self.network, mapping=relabel_nodes)

    @classmethod
    def name(cls):
        return cls.__name__

    @abstractmethod
    def load_network(self, file_resources: Dict, source_col_name: str, target_col_name: str, edge_attr: List[str],
                     directed: bool, filters: Dict[str, Any], blocksize=None):
        """
        Handles data processing from `file_resources` to a Pandas DataFrame which contain edgelist data, then constructs
        and return a NetworkX Graph.
        Args:
            file_resources: a dict of file name and file path/object
            source_col_name (str): column name of the dataframe for source in the edge
            target_col_name (str): column name of the dataframe for target in the edge
            edge_attr (list): list of str for column data to include in each edge
            directed (bool): True to return a DiGraph(), else Graph()
            filters: A dict of {column name: column values} to filter the dataframe
            blocksize ():
        Returns:
            network: a NetworkX Graph or DiGraph
        """
        raise NotImplementedError

    def get_interactions(self, nodelist=None, data=False, inclusive=True, relabel_nodes: Dict[str, str] = None):
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
        if not hasattr(self, "network"):
            raise Exception(
                "{} does not have network interaction data yet. Must run load_network() and assign self.network field first.".format(
                    self.name()))

        g = self.network
        if relabel_nodes:
            g = nx.relabel_nodes(g, relabel_nodes)

        if nodelist is None:
            return g.edges(data=data)

        if inclusive:
            return g.subgraph(nodelist).edges(data=data)
        else:
            return g.edges(nbunch=nodelist, data=data)


class STRING(Interactions, SequenceDatabase):
    """Loads the STRING database from https://string-db.org/ .

    Default path: "https://stringdb-static.org/download/" .
    Default file_resources: {
        "protein.actions.txt": f"protein.actions.v11.0/{species_id}.protein.actions.v11.0.txt.gz",
        "protein.links.txt": f"protein.links.v11.0/{species_id}.protein.links.v11.0.txt.gz",
        "protein.info.txt": f"protein.info.v11.0/{species_id}.protein.info.v11.0.txt.gz",
        "protein.sequences.fa": f"protein.sequences.v11.0/{species_id}.protein.sequences.v11.0.fa.gz"
    }

    Edge attributes for protein.actions.txt include ["mode", 'action', 'is_directional', 'a_is_acting' "score"]
    Edge attributes for protein.actions.txt include ["combined_score"]
    """
    COLUMNS_RENAME_DICT = {
        "#string_protein_id": "protein_id",
        "protein_external_id": "protein_id",
        "preferred_name": "protein_name",
        'combined_score': "weight",
    }

    def __init__(self, path="https://stringdb-static.org/download/", file_resources=None,
                 species_id: Union[str, List[str]] = "9606", version="v11.0",
                 source_col_name="protein1", target_col_name="protein2",
                 source_index="protein_name", target_index="protein_name",
                 edge_attr=['weight'], directed=False,
                 relabel_nodes=None, index_col='#string_protein_id', blocksize=None, **kwargs):
        """

        Args:
            path ():
            file_resources ():
            species_id (): List of str of species id's
                Must provide species_id number to download the species-specific STRING dataset.

            version ():
            source_col_name ():
            target_col_name ():
            source_index ():
            target_index ():
            edge_attr ():
            directed ():
            relabel_nodes ():
            verbose ():
            blocksize ():
        """
        self.version = version
        self.species_id = copy.copy(species_id)

        if file_resources is None:
            file_resources = {}
            if isinstance(species_id, (Iterable, str)) and len(species_id):
                species_list = species_id if isinstance(species_id, list) else [species_id]
                for id in species_list:
                    file_resources[f"{id}.protein.info.txt.gz"] = \
                        os.path.join(path, f"protein.info.{version}/{id}.protein.info.{version}.txt.gz")
                    file_resources[f"{id}.protein.links.txt.gz"] = \
                        os.path.join(path, f"protein.links.{version}/{id}.protein.links.{version}.txt.gz")
                    file_resources[f"{id}.protein.sequences.fa.gz"] = \
                        os.path.join(path, f"protein.sequences.{version}/{id}.protein.sequences.{version}.fa.gz")
            else:
                file_resources["protein.info.txt.gz"] = os.path.join(path, f"protein.info.{version}.txt.gz")
                file_resources["protein.links.txt.gz"] = os.path.join(path, f"protein.links.{version}.txt.gz")
                file_resources["protein.sequences.fa.gz"] = os.path.join(path, f"protein.sequences.{version}.fa.gz")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, directed=directed, relabel_nodes=relabel_nodes,
                         index_col=index_col,
                         blocksize=blocksize, **kwargs)

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        # Load nodes
        data_dfs = []
        if blocksize:
            for filename in [fn for fn, path in file_resources.items() \
                             if 'info.txt' in fn and isinstance(path, str)]:
                compression = 'gzip' if filename.endswith(".gz") else None
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    df = dd.read_table(file_resources[filename], na_values=['annotation not available'],
                                       low_memory=True, compression=compression,
                                       dtype={'protein_size': 'int8'},
                                       blocksize=None if isinstance(blocksize, bool) else blocksize)

                if self.keys is not None:
                    df = df.loc[df[self.index_col].isin(self.keys)]

                if self.index_col:
                    df = df.set_index(self.index_col, sorted=True)
                data_dfs.append(df)
        else:
            for filename in file_resources:
                if filename.endswith("protein.info.txt"):
                    df = pd.read_table(file_resources[filename], na_values=['annotation not available'],
                                       dtype={'protein_size': 'int8'},
                                       index_col=self.index_col, low_memory=True)
                    data_dfs.append(df)

        # Need to load .data df first
        nodes_df = dd.concat(data_dfs, axis=0, interleave_partitions=True) \
            if blocksize else pd.concat(data_dfs, axis=0)
        logger.info(f"nodes_df {nodes_df.columns}")
        nodes_df = nodes_df.rename(columns=self.COLUMNS_RENAME_DICT)
        if nodes_df.index.name in self.COLUMNS_RENAME_DICT:
            nodes_df.index = nodes_df.index.rename(self.COLUMNS_RENAME_DICT[nodes_df.index.name])

        return nodes_df

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        # Load edges
        edges_dfs = []
        if blocksize:
            for filename in [fn for fn, path in file_resources.items() \
                             if "links.txt" in fn and isinstance(path, str)]:
                compression = 'gzip' if filename.endswith(".gz") else None
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    df: dd.DataFrame = dd.read_table(file_resources[filename], sep=" ", low_memory=True,
                                                     compression=compression,
                                                     blocksize=None if isinstance(blocksize, bool) else blocksize)
                if compression:
                    logger.info(f"Repartitioning {filename} from {df.npartitions} "
                                f"partitions to {blocksize}-size partitions")
                    df = df.repartition(partition_size=blocksize)

                if self.keys is not None:
                    df = df.loc[df['protein1'].isin(self.keys) | df['protein2'].isin(self.keys)]
                edges_dfs.append(df)
        else:
            for filename in file_resources:
                if filename.endswith("links.txt"):
                    df = pd.read_table(file_resources[filename], sep=" ", low_memory=True)
                    edges_dfs.append(df)

        edges_df = dd.concat(edges_dfs, axis=0) if blocksize else pd.concat(edges_dfs, axis=0)
        edges_df = edges_df.rename(columns=self.COLUMNS_RENAME_DICT)
        logger.info(f"{self.name()}-{self.species_id}: {edges_df.columns.tolist()}, {edges_df.shape}")

        # Determine which edge attr to add
        if isinstance(edge_attr, (list, tuple)):
            cols = edges_df.columns.intersection(edge_attr + [source_col_name, target_col_name])
            edges_df = edges_df[cols]
            use_attrs = True
        else:
            use_attrs = False

        # return None
        # Add edgelist to nx Graph
        if blocksize:
            def edgelist2graph(edgelist_df: pd.DataFrame) -> nx.Graph():
                if edgelist_df.shape[0] == 1 and edgelist_df.iloc[0, 0] == 'foo': return
                G = nx.from_pandas_edgelist(edgelist_df, source=source_col_name, target=target_col_name,
                                            edge_attr=use_attrs, create_using=nx.DiGraph() if directed else nx.Graph())
                return G

            # Add edges to Networkx Graph
            edges_df: dd.DataFrame
            graphs_df = edges_df.map_partitions(edgelist2graph, meta=pd.Series([nx.Graph])).compute()
            G = nx.compose_all(graphs_df.values)

        else:
            G = nx.from_pandas_edgelist(edges_df, source=source_col_name, target=target_col_name,
                                        edge_attr=use_attrs, create_using=nx.DiGraph() if directed else nx.Graph())

        # self.protein_id2name = nodes_df.set_index("protein_id")["protein_name"].to_dict()
        # network = nx.relabel_nodes(network, self.protein_id2name)
        return G

    def get_sequences(self, index="protein_id", omic=None, agg=None):
        if hasattr(self, "seq_dict"):
            return self.seq_dict

        self.seq_dict = {}
        collisions = 0
        for record in SeqIO.parse(self.file_resources["protein.sequences.fa"], "fasta"):
            gene_id = str(record.name)

            sequence_str = str(record.seq)
            if index == "protein_name":
                key = self.protein_id2name[gene_id]
            elif index == "protein_id":
                key = gene_id

            if key in self.seq_dict:
                collisions += 1

            self.seq_dict[key] = sequence_str

        logger.warn("Seq {} collisions: {}".format(index, collisions))
        return self.seq_dict


class GeneMania(Interactions):
    """Loads the GeneMania database from  .

    Default path: local_directory .
    Default file_resources: {
        "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt": "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt",
        "identifier_mappings.txt": "identifier_mappings.txt",
    }
    """

    def __init__(self, path, file_resources=None, source_col_name="Gene_A", target_col_name="Gene_B",
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=None, filters=None, directed=True, relabel_nodes=None, **kwargs):
        if edge_attr is None:
            edge_attr = ["Weight"]
        if file_resources is None:
            file_resources = {}
            file_resources["COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt"] = os.path.join(path,
                                                                                        "COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt")
            file_resources["identifier_mappings.txt"] = os.path.join(path,
                                                                     "identifier_mappings.txt")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, filters=filters, directed=directed, relabel_nodes=relabel_nodes, **kwargs)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        interactions = pd.read_table(file_resources["COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt"], low_memory=True)
        identifier = pd.read_table(file_resources["identifier_mappings.txt"])

        # Rename ENSG ID's to gene names
        identifier = identifier[identifier["Source"] == "Gene Name"]
        id_mapping = pd.Series(identifier["Name"].values, index=identifier["Preferred_Name"]).to_dict()
        interactions.replace(id_mapping, inplace=True)

        genemania_RNA_RNA_network = nx.from_pandas_edgelist(interactions, source=source_col_name,
                                                            target=target_col_name,
                                                            edge_attr=edge_attr,
                                                            create_using=nx.DiGraph())
        return genemania_RNA_RNA_network


class IntAct(Interactions):

    def __init__(self, path, file_resources: Dict, source_col_name: str = None, target_col_name: str = None,
                 source_index: str = None, target_index: str = None, edge_attr: List[str] = None, filters: dict = None,
                 directed: bool = True, relabel_nodes: dict = None, blocksize=None, **kwargs):
        super().__init__(path, file_resources, source_col_name, target_col_name, source_index, target_index, edge_attr,
                         filters, directed, relabel_nodes, blocksize, **kwargs)


class BioGRID(Interactions):
    """Loads the BioGRID database from https://thebiogrid.org .

    Default path: "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/" .
    Default file_resources: {
        "BIOGRID-ALL-LATEST.tab2.zip": "BIOGRID-ALL-LATEST.tab2.zip",
    }
    """

    def __init__(self, path="https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/",
                 file_resources=None, source_col_name="Official Symbol Interactor A",
                 target_col_name="Official Symbol Interactor B",
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=['Score', 'Throughput', 'Experimental System', 'Experimental System Type'],
                 filters=None, directed=False, relabel_nodes=None, **kwargs):
        """

        Args:
            path ():
            file_resources ():
            source_col_name ():
            target_col_name ():
            source_index ():
            target_index ():
            edge_attr ():
            filters (): Default None, example {"Organism Interactor A": 9606}.
            directed ():
            relabel_nodes ():
            **kwargs ():
        """
        if file_resources is None:
            file_resources = {}
            file_resources["BIOGRID-ALL-LATEST.tab2.zip"] = os.path.join(path, "BIOGRID-ALL-LATEST.tab2.zip")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, filters=filters, directed=directed, relabel_nodes=relabel_nodes, **kwargs)

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        args = dict(na_values=["-"], header=0, low_memory=True,
                    # usecols=['Official Symbol Interactor A', 'Official Symbol Interactor B',
                    #          'Organism Interactor A', 'Score', 'Throughput', 'Qualifications',
                    #          'Modification', 'Phenotypes', 'Source Database'],
                    dtype={'Score': 'float', 'Entrez Gene Interactor A': 'str', 'Entrez Gene Interactor B': 'str',
                           'Organism Interactor A': 'category', 'Organism Interactor B': 'category'})

        if blocksize:
            df = dd.read_table(file_resources["BIOGRID-ALL-LATEST.tab2"], blocksize=blocksize, **args, )
        else:
            df = pd.read_table(file_resources["BIOGRID-ALL-LATEST.tab2"], **args, )

        return df

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        return
        df = self.data
        df = filter_rows(df, filters)
        network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                          edge_attr=edge_attr,
                                          create_using=nx.DiGraph() if directed else nx.Graph())
        return network


class MiRTarBase(Interactions):
    """Loads the  database from  .

            Default path:  .
            Default file_resources: {
                "": "",
                "": "",
                "": "",
            }
            """

    def __init__(self, path="http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/", file_resources=None,
                 source_col_name="miRNA", target_col_name="Target Gene",
                 source_index="transcript_name", target_index="gene_name",
                 edge_attr=None,
                 filters=None,
                 directed=True,
                 relabel_nodes=None,
                 strip_mirna_name=False, **kwargs):
        """

        Args:
            path ():
            file_resources ():
            source_col_name ():
            target_col_name ():
            source_index ():
            target_index ():
            edge_attr ():
            filters (): default None, example {"Species (Target Gene)": "Homo sapiens"}
            directed ():
            relabel_nodes ():
            strip_mirna_name ():
            **kwargs ():
        """
        if edge_attr is None:
            edge_attr = ["Support Type"]
        self.strip_mirna_name = strip_mirna_name

        if file_resources is None:
            file_resources = {}
            file_resources["miRTarBase_MTI.xlsx"] = os.path.join(path, "miRTarBase_MTI.xlsx")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, filters=filters, directed=directed, relabel_nodes=relabel_nodes, **kwargs)

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        df = pd.read_excel(self.file_resources["miRTarBase_MTI.xlsx"])
        return df

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        df = self.data
        df = filter_rows(df, filters)

        if self.strip_mirna_name:
            df['miRNA'] = df['miRNA'].str.lower()
            df['miRNA'] = df['miRNA'].str.replace("-3p.*|-5p.*", "", regex=True)

        mir_target_network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                                     edge_attr=edge_attr,
                                                     create_using=nx.DiGraph() if directed else nx.Graph())
        return mir_target_network


class LncBase(Interactions, Database):
    """Loads the LncBase database from http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=lncbasev2%2Findex .

    Default path: local_directory .
    Default file_resources: {
        "LncBasev2_download.csv": "LncBasev2_download.csv"",
    }
    """

    def __init__(self, path, file_resources=None, strip_mirna_name=False,
                 source_col_name="mirna", target_col_name="geneId",
                 source_index="transcript_name", target_index="gene_id",
                 edge_attr=None,
                 filters=None,
                 directed=True,
                 relabel_nodes=None, ):
        """

        Args:
            path ():
            file_resources ():
            strip_mirna_name ():
            source_col_name ():
            target_col_name ():
            source_index ():
            target_index ():
            edge_attr ():
            filters (): default None. Example: {"species": "Homo sapiens"}
            directed ():
            relabel_nodes ():
        """
        self.strip_mirna_name = strip_mirna_name

        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        if file_resources is None:
            file_resources = {}
            file_resources["LncBasev2_download.csv"] = os.path.join(path, "LncBasev2_download.csv")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name,
                         source_index=source_index, target_index=target_index, edge_attr=edge_attr, filters=filters,
                         directed=directed, relabel_nodes=relabel_nodes)

    def get_rename_dict(self, from_index="geneId", to_index="geneName"):
        lncbase_df = pd.read_table(self.file_resources["LncBasev2_download.csv"], low_memory=True)
        gene_id_to_gene_name_dict = pd.Series(lncbase_df["geneName"].values,
                                              index=lncbase_df["geneId"]).to_dict()
        return gene_id_to_gene_name_dict

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        df = pd.read_table(file_resources["LncBasev2_download.csv"], low_memory=True)
        df.replace({"species": {"Homo Sapiens": "Homo sapiens", "Mus Musculus": "Mus musculus"}}, inplace=True)
        return df

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        df = self.data
        df = filter_rows(df, filters)

        if self.strip_mirna_name:
            df['mirna'] = df['mirna'].str.lower()
            df['mirna'] = df['mirna'].str.replace("-3p.*|-5p.*", "", regex=True)

        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                                               edge_attr=edge_attr,
                                                               create_using=nx.DiGraph() if directed else nx.Graph())
        return lncBase_lncRNA_miRNA_network


class TargetScan(Interactions, Database):
    """Loads the TargetScan database from "http://www.targetscan.org/" .

    Default path: "http://www.targetscan.org/vert_72/vert_72_data_download/" .
    Default file_resources: {
        "miR_Family_Info.txt": "miR_Family_Info.txt.zip",
        "Predicted_Targets_Info.default_predictions.txt": "Predicted_Targets_Info.default_predictions.txt.zip",
        "": "",
    }
    """

    def __init__(self, path="http://www.targetscan.org/vert_72/vert_72_data_download/", file_resources=None,
                 source_col_name="MiRBase ID", target_col_name="Gene Symbol",
                 source_index="transcript_name", target_index="transcript_name",
                 edge_attr=["tissue", "positive_negative"], directed=True, relabel_nodes=None, species=9606,
                 strip_mirna_name=False):
        if edge_attr is None:
            edge_attr = ["tissue", "positive_negative"]
        self.strip_mirna_name = strip_mirna_name
        self.species = species
        if file_resources is None:
            file_resources = {}
            file_resources["miR_Family_Info.txt.zip"] = os.path.join(path, "miR_Family_Info.txt.zip")
            file_resources["Predicted_Targets_Info.default_predictions.txt"] = os.path.join(path,
                                                                                            "Predicted_Targets_Info.default_predictions.txt")

        super().__init__(path=path, file_resources=file_resources)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
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

        family_to_miR_df.set_genes_index("miR Family", inplace=True)
        family_interactions_df.set_genes_index("miR Family", inplace=True)
        mir_interactions_df = family_interactions_df.join(family_to_miR_df, how='outer', on="miR Family").reset_index()

        # Standardize MiRBase ID to miRNA names obtained from RNA-seq hg19
        if self.strip_mirna_name:
            mir_interactions_df['MiRBase ID'] = mir_interactions_df['MiRBase ID'].str.lower()
            mir_interactions_df['MiRBase ID'] = mir_interactions_df['MiRBase ID'].str.replace("-3p.*|-5p.*", "")

        return mir_interactions_df


class LncReg(Interactions):
    """Loads the  database from  .

    Default path:  .
    Default file_resources: {
        "": "",
        "": "",
        "": "",
    }
    """
    def __init__(self, path, file_resources,
                 source_col_name='A_name_in_paper', target_col_name='B_name_in_paper',
                 source_index="transcript_name", target_index="gene_name",
                 edge_attr=["relationship", "mechanism", "pmid"], filters=None, directed=True, relabel_nodes=None,
                 verbose=False):
        if file_resources is None:
            file_resources = {}
            file_resources["data.xlsx"] = os.path.join(path, "data.xlsx")

        super().__init__(path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, filters=filters, directed=directed, relabel_nodes=relabel_nodes,
                         verbose=verbose)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        df = pd.read_excel(self.file_resources["data.xlsx"])
        print(self.name(), df.columns.tolist())

        df = df[df["species"] == "Homo sapiens"]
        df.loc[df["B_category"] == "miRNA", "B_name_in_paper"] = df[df["B_category"] == "miRNA"][
            "B_name_in_paper"].str.replace("-3p.*|-5p.*", "")
        df.loc[df["B_category"] == "miRNA", "B_name_in_paper"] = df[df["B_category"] == "miRNA"][
            "B_name_in_paper"].str.replace("MIR", "hsa-mir-")
        df.loc[df["B_category"] == "miRNA", "B_name_in_paper"] = df[df["B_category"] == "miRNA"][
            "B_name_in_paper"].str.replace("let-", "hsa-let-")

        LncReg_lncRNA_RNA_network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                                            edge_attr=edge_attr,
                                                            create_using=nx.DiGraph())
        return LncReg_lncRNA_RNA_network


class lncRInter(Interactions):
    """Loads the  database from  .

    Default path:  .
    Default file_resources: {
        "": "",
        "": "",
        "": "",
    }
    """
    def __init__(self, path, file_resources=None, source_col_name="lncrna",
                 target_col_name='Interacting partner',
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=None, filters=None,
                 directed=True, relabel_nodes=None):
        if edge_attr is None:
            edge_attr = ["Interaction Class", "Interaction Mode", "Tissue", "Phenotype"]
        if file_resources is None:
            file_resources = {}
            file_resources["human_interactions.txt"] = os.path.join(path, "human_interactions.txt")

        super().__init__(path, file_resources, source_col_name=source_col_name, target_col_name=target_col_name,
                         source_index=source_index, target_index=target_index, edge_attr=edge_attr, filters=filters,
                         directed=directed, relabel_nodes=relabel_nodes)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        lncRInter_df = pd.read_table(file_resources["human_interactions.txt"])
        print(self.name(), lncRInter_df.columns.tolist())

        lncRInter_df = filter_rows(lncRInter_df, filters)
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
    """Loads the  database from  .

            Default path:  .
            Default file_resources: {
                "": "",
                "": "",
                "": "",
            }
            """

    def __init__(self, path="http://123.59.132.21/lncrna2target/data/", file_resources=None, source_index="gene_name",
                 target_index="gene_name", edge_attr=None,
                 filters=None,
                 directed=True, relabel_nodes=None, version="high_throughput", ):
        """

        Args:
            filters (): default None, example {"species_id": 9606, "Species": "Homo sapiens"}.
            version (str): one of ["high_throughput", "low_throughput"].
                The high_throughput version of lncRNA2Target database is v2.0 and low_throughput is v1.0, according to the database's website.
            species_id (str, int): one of [9606, "Homo sapiens"].
                The species column in high_throughput is formatted in int (e.g. 9606) and in low_throughput is in str (e.g. "Homo sapiens")
        """
        self.version = version
        if file_resources is None:
            file_resources = {}
            file_resources["lncRNA_target_from_high_throughput_experiments.txt.rar"] = os.path.join(path,
                                                                                                    "lncrna_target.rar")
            file_resources["lncRNA_target_from_low_throughput_experiments.xlsx"] = os.path.join(path,
                                                                                                "lncRNA_target_from_low_throughput_experiments.xlsx")

        if self.version == "high_throughput":
            super().__init__(path, file_resources, source_col_name="lncrna_symbol", target_col_name="gene_symbol",
                             source_index=source_index, target_index=target_index, edge_attr=edge_attr, filters=filters,
                             directed=directed, relabel_nodes=relabel_nodes)
        if self.version == "low_throughput":
            super().__init__(path, file_resources, source_col_name="GENCODE_gene_name",
                             target_col_name="Target_official_symbol", source_index=source_index,
                             target_index=target_index, edge_attr=edge_attr, filters=filters, directed=directed,
                             relabel_nodes=relabel_nodes)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        if self.version == "high_throughput":
            return self.load_network_high_throughput(file_resources, source_col_name, target_col_name, edge_attr,
                                                     directed)
        elif self.version == "low_throughput":
            return self.load_network_low_throughput(file_resources, source_col_name, target_col_name, edge_attr,
                                                    directed)
        else:
            raise Exception("LncRNA2Target version argument must be one of 'high_throughput' or 'low_throughput'")

    def load_network_high_throughput(self, file_resources, source_col_name="lncrna_symbol",
                                     target_col_name="gene_symbol",
                                     edge_attr=None, directed=True, filters=None):
        table = pd.read_table(file_resources["lncRNA_target_from_high_throughput_experiments.txt"], sep="\t")
        table = filter_rows(table, filters)
        print(self.name(), table.columns.tolist())

        table["lncrna_symbol"] = table["lncrna_symbol"].str.upper()
        table["lncrna_symbol"] = table["lncrna_symbol"].str.replace("LINC", "")
        table["gene_symbol"] = table["gene_symbol"].str.upper()
        lncrna2target_high_throughput_network = nx.from_pandas_edgelist(table,
                                                                        source=source_col_name,
                                                                        target=target_col_name,
                                                                        edge_attr=edge_attr,
                                                                        create_using=nx.DiGraph() if directed else nx.Graph())
        return lncrna2target_high_throughput_network

    def load_network_low_throughput(self, file_resources, source_col_name="GENCODE_gene_name",
                                    target_col_name="Target_official_symbol",
                                    edge_attr=None, directed=True, filters=None):
        table = pd.read_excel(file_resources["lncRNA_target_from_low_throughput_experiments.xlsx"])
        table = filter_rows(table, filters)
        print(self.name(), table.columns.tolist())

        table["Target_official_symbol"] = table["Target_official_symbol"].str.replace("(?i)(mir)", "hsa-mir-",
                                                                                      regex=True)
        table["Target_official_symbol"] = table["Target_official_symbol"].str.replace("--", "-")
        table["Target_official_symbol"].apply(lambda x: x.lower() if "mir" in x.lower() else x.upper())
        table["GENCODE_gene_name"] = table["GENCODE_gene_name"].str.upper()
        lncrna2target_low_throughput_network = nx.from_pandas_edgelist(table,
                                                                       source=source_col_name,
                                                                       target=target_col_name,
                                                                       edge_attr=edge_attr,
                                                                       create_using=nx.DiGraph() if directed else nx.Graph())
        return lncrna2target_low_throughput_network


class lncRNome(Interactions, Database):
    """Loads the lncRNome database from  .

    Default path:  .
    Default file_resources: {
        "": "",
        "": "",
        "": "",
    }
    """

    def __init__(self, path, file_resources, source_col_name='Gene Name', target_col_name='Binding miRNAs',
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=["miRNA Interaction Site", "Transcript ID"], directed=True, relabel_nodes=None,
                 blocksize=0):
        if file_resources is None:
            file_resources = {}
            file_resources["miRNA_binding_sites.txt"] = os.path.join(path, "miRNA_binding_sites.txt")
            file_resources["general_information.txt"] = os.path.join(path, "general_information.txt")

        super().__init__(path, file_resources=file_resources)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        df = pd.read_table(self.file_resources["miRNA_binding_sites.txt"], header=0)
        print(self.name(), df.columns.tolist())

        df['Binding miRNAs'] = df['Binding miRNAs'].str.lower()
        df['Binding miRNAs'] = df['Binding miRNAs'].str.replace("-3p.*|-5p.*", "", regex=True)

        lncRNome_miRNA_binding_sites_network = nx.from_pandas_edgelist(df, source=source_col_name,
                                                                       target=target_col_name,
                                                                       edge_attr=edge_attr,
                                                                       create_using=nx.DiGraph())

        return lncRNome_miRNA_binding_sites_network

    def load_dataframe(self, file_resources, blocksize=None):
        return pd.read_table(self.file_resources["general_information.txt"], header=0,
                             usecols=["Gene Name", "Transcript Name", "Transcript Type", "Location", "Strand"])


class NPInter(Interactions):
    """Loads the NPInter database from http://bigdata.ibp.ac.cn/npinter4/ .

    Default path: "http://bigdata.ibp.ac.cn/npinter4/download/" .
    Default file_resources: {
        "interaction_NPInterv4.expr.txt": "file/interaction_NPInterv4.expr.txt.gz",
    }
    """
    def __init__(self, path="http://bigdata.ibp.ac.cn/npinter4/download/", file_resources=None,
                 source_col_name='ncName', target_col_name='tarName',
                 source_index="gene_name", target_index="gene_name",
                 edge_attr=["tarType", "tissueOrCell", "tag", 'class', "level"],
                 filters=None,
                 directed=True, relabel_nodes=None, verbose=False):
        if file_resources is None:
            file_resources = {}
            file_resources["interaction_NPInterv4.expr.txt.gz"] = \
                os.path.join(path, "file/interaction_NPInterv4.expr.txt.gz")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, source_index=source_index, target_index=target_index,
                         edge_attr=edge_attr, filters=filters, directed=directed, relabel_nodes=relabel_nodes,
                         verbose=verbose)

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        df = pd.read_table(file_resources["interaction_NPInterv4.expr.txt"], header=0, na_values=["-"])
        print(self.name(), df.columns.tolist())
        df["ncName"] = df["ncName"].str.upper()
        df["ncName"] = df["ncName"].str.strip("LNCRNA-")
        df["ncName"] = df["ncName"].str.replace("MALAT-1", "MALAT1")
        df["ncName"] = df["ncName"].str.replace("^MIR-", "hsa-mir-", regex=True)
        df["ncName"] = df["ncName"].str.replace("^MICRORNA-", "hsa-mir-", regex=True)

        df["tarName"] = df["tarName"].str.upper()

        return df

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        df = self.data
        df = filter_rows(df, filters)

        lncRNome_miRNA_binding_sites_network = nx.from_pandas_edgelist(df, source=source_col_name,
                                                                       target=target_col_name,
                                                                       edge_attr=edge_attr,
                                                                       create_using=nx.DiGraph() if directed else nx.Graph())

        return lncRNome_miRNA_binding_sites_network


class StarBase(Interactions):
    """Loads the  database from  .

    Default path:  .
    Default file_resources: {
        "": "",
        "": "",
        "": "",
    }
    """

    def __init__(self, path, file_resources, source_col_name="geneName", target_col_name="pairGeneName",
                 source_index="gene_name", target_index="gene_name",
                 min_interactionNum=1, min_expNum=1,
                 edge_attr=None, directed=True, relabel_nodes=None, blocksize=0):
        if file_resources is None:
            file_resources = {}
            file_resources["starbase_3.0_lncrna_rna_interactions.csv"] = \
                os.path.join(path, "starbase_3.0_lncrna_rna_interactions.csv")
        self.min_interactionNum = min_interactionNum
        self.min_expNum = min_expNum
        super().__init__(path, file_resources, source_col_name, target_col_name, source_index, target_index, edge_attr,
                         directed, relabel_nodes, blocksize)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        df = pd.read_csv(self.file_resources["starbase_3.0_lncrna_rna_interactions.csv"], header=0)

        df.loc[df["pairGeneType"] == "miRNA", "pairGeneName"] = df[df["pairGeneType"] == "miRNA"][
            "pairGeneName"].str.lower()
        df.loc[df["pairGeneType"] == "miRNA", "pairGeneName"] = df[df["pairGeneType"] == "miRNA"][
            "pairGeneName"].str.replace("-3p.*|-5p.*", "")
        df = df[df["interactionNum"] >= self.min_interactionNum]
        df = df[df["expNum"] >= self.min_expNum]

        self.starBase_RNA_RNA_network = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                                                edge_attr=["interactionNum"],
                                                                create_using=nx.DiGraph())
        return self.starBase_RNA_RNA_network

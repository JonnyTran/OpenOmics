import copy
import gc
import os
from abc import abstractmethod
from collections.abc import Iterable
from typing import List, Dict, Any, Union, Optional

import dask.dataframe as dd
import networkx as nx
import pandas as pd
import scipy.sparse as ssp
from Bio import SeqIO
from logzero import logger

from openomics.database.base import Database
from openomics.database.sequence import SequenceDatabase, UniProt
from openomics.transforms.df import filter_rows

__all__ = ['STRING', 'GeneMania', 'IntAct', 'BioGRID', 'MiRTarBase', 'LncBase', 'TargetScan', 'TarBase',
           'LncReg', 'LncRNA2Target', 'lncRNome', 'NPInter', 'RNAInter', 'StarBase']

class Interactions(Database):
    edges: Optional[Union[pd.DataFrame, dd.DataFrame]]
    def __init__(self, path, file_resources: Dict, source_col_name: str = None, target_col_name: str = None,
                 edge_attr: List[str] = None, filters: Union[str, Dict[str, Union[str, List[str]]]] = None,
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
        self.filters = filters

        super().__init__(path=path, file_resources=file_resources, blocksize=blocksize, **kwargs)
        self.network = self.load_network(file_resources=self.file_resources, source_col_name=source_col_name,
                                         target_col_name=target_col_name, edge_attr=edge_attr, directed=directed,
                                         filters=filters, blocksize=blocksize)

        if relabel_nodes is not None:
            self.network = nx.relabel_nodes(self.network, mapping=relabel_nodes)

        self.close()

    @classmethod
    def name(cls):
        return cls.__name__

    @abstractmethod
    def load_network(self, file_resources: Dict, source_col_name: str, target_col_name: str,
                     edge_attr: Union[str, List[str]], directed: bool, filters: Dict[str, Any], blocksize=None) \
        -> nx.Graph:
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
            g = nx.relabel_nodes(g, relabel_nodes, copy=False)

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
        "{species_id}.protein.info.txt.gz": f"protein.info.{version}/{species_id}.protein.info.{version}.txt.gz",
        "{species_id}.protein.aliases.txt.gz": f"protein.links.{version}/{species_id}.protein.aliases.{version}.txt.gz",
        "{species_id}.protein.links.txt.gz": f"protein.links.{version}/{species_id}.protein.links.{version}.txt.gz",
        "{species_id}.protein.sequences.fa.gz": f"protein.sequences.{version}/{species_id}.protein.sequences.{version}.fa.gz"
    }

    Edge attributes for protein.actions.txt include ["mode", 'action', 'is_directional', 'a_is_acting' "score"]
    Edge attributes for protein.actions.txt include ["combined_score"]
    """
    COLUMNS_RENAME_DICT = {
        "#string_protein_id": "string_protein_id",
        "protein_external_id": "protein_id",
        "preferred_name": "gene_name",
        '#ncbi_taxid': 'species_id',
        'string_protein_id_2': 'homologous_protein_id',
    }

    def __init__(self, path="https://stringdb-static.org/download/", file_resources=None,
                 species_id: Union[str, List[str]] = "9606", version="v11.0",
                 source_col_name="protein1", target_col_name="protein2",
                 edge_attr='combined_score', directed=False,
                 relabel_nodes=None,
                 index_col='#string_protein_id',
                 keys=None,
                 alias_types={'Ensembl_UniProt', 'Ensembl_UniProt_AC'},
                 blocksize=None, **kwargs):
        """

        Args:
            path ():
            file_resources ():
            species_id (): List of str of species id's
                Provide a species_id string or a list of species_id's to download the species-specific STRING dataset, and
                integrate them. If species_id is None, then download the full-dataset version of STRING, which is very
                time-consuming.
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
        self.alias_types = alias_types
        assert isinstance(edge_attr, str)

        if file_resources is None:
            file_resources = {}
            if isinstance(species_id, (Iterable, str)) and len(species_id):
                species_list = [species_id] if isinstance(species_id, str) else species_id
                for species in species_list:
                    file_resources[f"{species}.protein.info.txt.gz"] = \
                        os.path.join(path, f"protein.info.{version}/{species}.protein.info.{version}.txt.gz")
                    file_resources[f"{species}.protein.links.txt.gz"] = \
                        os.path.join(path, f"protein.links.{version}/{species}.protein.links.{version}.txt.gz")
                    file_resources[f"{species}.protein.links.detailed.txt.gz"] = \
                        os.path.join(path, f"protein.links.detailed.{version}/"
                                           f"{species}.protein.links.detailed.{version}.txt.gz")
                    file_resources[f"{species}.protein.homology.txt.gz"] = \
                        os.path.join(path, f"protein.homology.{version}/{species}.protein.homology.{version}.txt.gz")
                    file_resources[f"{species}.clusters.proteins.txt.gz"] = \
                        os.path.join(path, f"clusters.proteins.{version}/{species}.clusters.proteins.{version}.txt.gz")
                    file_resources[f"{species}.protein.aliases.txt.gz"] = \
                        os.path.join(path, f"protein.aliases.{version}/{species}.protein.aliases.{version}.txt.gz")
                    file_resources[f"{species}.enrichment.terms.txt.gz"] = \
                        os.path.join(path, f"enrichment.terms.{version}/{species}.enrichment.terms.{version}.txt.gz")
                    file_resources[f"{species}.protein.sequences.fa.gz"] = \
                        os.path.join(path, f"protein.sequences.{version}/{species}.protein.sequences.{version}.fa.gz")
            else:
                file_resources["protein.info.txt.gz"] = os.path.join(path, f"protein.info.{version}.txt.gz")
                file_resources["protein.links.txt.gz"] = os.path.join(path, f"protein.links.{version}.txt.gz")
                file_resources["protein.sequences.fa.gz"] = os.path.join(path, f"protein.sequences.{version}.fa.gz")
        else:
            if isinstance(self.species_id, Iterable):
                file_resources = {fn: fp for fn, fp in file_resources.items() \
                                  if any(fn.startswith(species) for species in self.species_id)}

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, edge_attr=edge_attr, directed=directed,
                         relabel_nodes=relabel_nodes, blocksize=blocksize, index_col=index_col, keys=keys,
                         col_rename=STRING.COLUMNS_RENAME_DICT, **kwargs)

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        # Load nodes
        dfs = []
        if blocksize:
            for filename in [fn for fn, path in file_resources.items() \
                             if 'info.txt' in fn and isinstance(path, str)]:
                compression = 'gzip' if filename.endswith(".gz") else None
                info_df = dd.read_table(file_resources[filename], na_values=['annotation not available'],
                                        low_memory=True, compression=compression,
                                        dtype={'protein_size': 'int8'},
                                        blocksize=None if isinstance(blocksize, bool) else blocksize)

                if self.keys is not None:
                    info_df = info_df.loc[info_df[self.index_col].isin(self.keys)]

                if self.index_col:
                    info_df = info_df.set_index(self.index_col, sorted=True)

                # Join other attributes to node_info
                species_id = filename.split(".")[0]
                attrs = self.load_accessory_data(file_resources, species_id=species_id,
                                                 alias_types=self.alias_types, blocksize=False)
                if attrs is not None:
                    new_cols = attrs.columns.difference(info_df.columns)
                    info_df = info_df.join(attrs[new_cols], on=self.index_col)

                dfs.append(info_df)
        else:
            for filename in file_resources:
                if filename.endswith("protein.info.txt"):
                    info_df = pd.read_table(file_resources[filename], na_values=['annotation not available'],
                                            dtype={'protein_size': 'int8'},
                                            index_col=self.index_col, low_memory=True)
                    index_split = info_df['#string_protein_id'].str.split(".", expand=True, n=1)
                    info_df = info_df.assign(species_id=index_split[0], protein_embl_id=index_split[1])

                    # Join other attributes to node_info
                    species_id = filename.split(".")[0]
                    attrs = self.load_accessory_data(file_resources, species_id=species_id,
                                                     alias_types=self.alias_types,
                                                     blocksize=blocksize)
                    if attrs is not None:
                        new_cols = attrs.columns.difference(info_df.columns)
                        info_df = info_df.join(attrs[new_cols], on=self.index_col)
                    dfs.append(info_df)

        if not len(dfs):
            raise Exception("Must provide at least one 'protein.info.txt' file.")

        if blocksize:
            protein_info: dd.DataFrame = dd.concat(dfs, axis=0, interleave_partitions=True)
        else:
            protein_info = pd.concat(dfs, axis=0)

        return protein_info

    def load_accessory_data(self, file_resources: Dict[str, str], species_id: str,
                            accessory_files=['protein.aliases', 'protein.homology', 'protein.enrichment',
                                             'clusters.proteins'],
                            alias_types={'Ensembl_UniProt', 'Ensembl_UniProt_AC'}, blocksize=False, ) \
        -> Union[pd.DataFrame, dd.DataFrame]:
        """
        Stack the annotations files for the provided `species_id`, such that rows in the annotations are filtered by
        `keys` (if not null), indexed by "#string_protein_id", and with attributes transformed to a dataframe columns.

        Args:
            file_resources (): a dict of filename and filepath
            species_id (str): the species_id string which is used to select only files that have the same prefix.
            accessory_files (List[str]):
                A list of strings that specify which types of annotation files to integrate, i.e., only select files
                having a substring matching one of these.
                Default ['protein.aliases', 'protein.homology', 'protein.enrichment', 'clusters.proteins'].
            alias_types (): a set of string, default {'Ensembl_UniProt_AC'}
                A set of `source` values in the `protein.aliases` annotation to aggregate `alias`'s for.
                Must be a subset of {'Ensembl_Source', 'Ensembl_gene', 'Ensembl_transcript', 'Ensembl_UniGene',
                    'Ensembl_RefSeq_short', 'Ensembl_RefSeq', 'Ensembl_OTTG', 'Ensembl_OTTP', 'Ensembl_UCSC',
                    'Ensembl_UniProt', 'Ensembl_UniProt_AC', 'Ensembl_EntrezGene', 'Ensembl_EMBL', 'Ensembl_protein_id'}
            blocksize (bool): Recommended to use Pandas to avoid uncessary overhead.

        Returns:
            dd.Dataframe or pd.DataFrame

        """
        allowed_prefixes = {'protein.aliases', 'protein.homology', 'protein.enrichment', 'clusters.proteins'}
        if not set(accessory_files).issubset(allowed_prefixes):
            logger.warn(f'{set(accessory_files).difference(allowed_prefixes)} files are not supported')

        select_files = []
        for fn, path in file_resources.items():
            if fn.startswith(species_id) and any(ftype in fn for ftype in accessory_files):
                select_files.append(fn)

        dfs = []
        for filename in select_files:
            args = dict(
                low_memory=True,
                dtype={'cluster_id': 'category', '#ncbi_taxid': 'category', 'category': 'category',
                       'source': 'category'})
            compression = 'gzip' if filename.endswith(".gz") else None
            if blocksize:
                if not isinstance(file_resources[filename], str): continue
                df = dd.read_table(file_resources[filename], compression=compression, **args)
            else:
                df = pd.read_table(file_resources[filename], **args)

            # Set index for df
            for col in ['#string_protein_id', 'protein_id', '#string_protein_1']:
                if col in df.columns:
                    df = df.set_index(col, sorted=True) if blocksize else df.set_index(col)
                    break

            # Set index
            if df.index.name is None:
                continue
            elif self.index_col and df.index.name != self.index_col:
                df.index = df.index.rename(self.index_col)
            if blocksize:
                assert df.known_divisions

            # Filter rows
            if self.keys is not None:
                df = df.loc[df.index.isin(self.keys)]

            # Groupby on index and perform appropriate transforms depending on the annotation type
            if 'protein.homology' in filename:
                df = df.loc[df.index != df['string_protein_id_2']]
                df = df.groupby(self.index_col)['string_protein_id_2'].unique().to_frame()
                # TODO ignored column of size of homologous regions

            elif 'clusters.protein' in filename:
                df = df.groupby(self.index_col)[['cluster_id', '#ncbi_taxid']].unique()

            elif 'protein.enrichment' in filename:
                df = df.groupby(self.index_col)['term'].unique().to_frame()

            elif 'protein.aliases' in filename:
                df = df.loc[df['source'].isin(alias_types)]
                df['source'] = df['source'].cat.set_categories(alias_types)
                if blocksize:
                    # Set alias values to lists so pivot_table(..., aggfunc='sum') will concatenate them
                    df = df.assign(alias=df['alias'].map(lambda x: [x], meta=pd.Series([[""]])))
                    df = dd.pivot_table(df.reset_index(),
                                        index='#string_protein_id', columns='source', values='alias', aggfunc='sum')
                else:
                    df = df.reset_index().groupby([self.index_col, 'source'])['alias'].unique().unstack(level=1)

            if blocksize and not df.known_divisions:
                df.divisions = df.compute_current_divisions()

            if not len(df.index):
                continue

            dfs.append(df)

        if dfs:
            attrs = dd.concat(dfs, axis=1) if blocksize else pd.concat(dfs, axis=1)
        else:
            attrs = None

        return attrs

    def load_network(self, file_resources, source_col_name='protein1', target_col_name='protein2',
                     edge_attr: str = 'combined_score', directed=False, filters=None, blocksize=None):
        keys = self.data.index.compute() if isinstance(self.data, dd.DataFrame) else self.data.index
        select_files = [fn for fn, path in file_resources.items() if "links" in fn]

        # Load edges
        edges_dfs = []
        for filename in select_files:
            args = dict(sep=" ", low_memory=True,
                        dtype={'protein1': 'category', 'protein2': 'category',
                               'neighborhood': 'uint8', 'fusion': 'uint8', 'cooccurence': 'uint8',
                               'coexpression': 'uint8', 'experimental': 'uint8', 'database': 'uint8',
                               'textmining': 'uint8', 'combined_score': 'uint8'})
            if blocksize:
                if not isinstance(file_resources[filename], str): continue
                compression = 'gzip' if filename.endswith(".gz") else None
                df: dd.DataFrame = dd.read_table(file_resources[filename], compression=compression, **args,
                                                 blocksize=None if isinstance(blocksize, bool) else blocksize)

                if compression:
                    logger.info(f"Repartitioning {filename} from {df.npartitions} "
                                f"partitions to {blocksize}-size partitions")
                    df = df.repartition(partition_size=blocksize)

            else:
                df = pd.read_table(file_resources[filename], **args)

            df = df.loc[df[source_col_name].isin(keys) & df[target_col_name].isin(keys)]
            edges_dfs.append(df)

        if len(edges_dfs) == 0:
            return

        # Concatenate multiple edgelists into dataframe
        edges_df = dd.concat(edges_dfs, axis=0) if blocksize else pd.concat(edges_dfs, axis=0)
        edges_df = edges_df.rename(columns=self.COLUMNS_RENAME_DICT)
        logger.info(f"{self.name()}-{self.species_id}: {edges_df.columns.tolist()}, {edges_df.shape}")

        # Convert edge_attr (edge weights) from 3 digit integer to float
        assignfunc = {}
        for col in (edge_attr if isinstance(edge_attr, Iterable) else [edge_attr]):
            if col in edges_df.columns:
                assignfunc[col] = edges_df[col].astype('float16') / 1000
        if assignfunc:
            edges_df = edges_df.assign(**assignfunc)

        edges_df = filter_rows(edges_df, filters=filters)

        self.edges = edges_df
        # Set ordering for rows and columns
        node2idx = {node: i for i, node in enumerate(keys)}

        if blocksize:
            edges_df: dd.DataFrame

            def edgelist2adj(df: pd.DataFrame) -> ssp.coo_matrix:
                if df.shape[0] == 1 and df.iloc[0, 0] == 'foo':
                    return None

                df = df.assign(row=df[source_col_name].map(node2idx).astype('int'),
                               col=df[target_col_name].map(node2idx).astype('int'))
                df = df.dropna(subset=['row', 'col'])

                if df.shape[0] == 0:
                    return None

                coo_adj = ssp.coo_matrix((df[edge_attr], (df['row'], df['col'])),
                                         shape=(len(keys), len(keys)))
                return coo_adj

            # Create a sparse adjacency matrix each partition, then combine them
            adj = edges_df.reduction(chunk=edgelist2adj,
                                     aggregate=lambda x: x.dropna().sum() if not x.isna().all() else None,
                                     meta=pd.Series([ssp.coo_matrix])).compute()
            assert len(adj) == 1, f"len(adj) = {len(adj)}"

            G = nx.from_scipy_sparse_matrix(adj[0], create_using=nx.DiGraph() if directed else nx.Graph(),
                                            edge_attribute=edge_attr)

            idx2node = {i: node for i, node in enumerate(keys)}
            G = nx.relabel_nodes(G, mapping=idx2node, copy=False)
            del adj
            gc.collect()

        else:
            # Determine which edge attr to add
            if isinstance(edge_attr, (list, tuple)):
                cols = edges_df.columns.intersection(edge_attr + [source_col_name, target_col_name])
                edges_df = edges_df[cols]
                use_attrs = True
            elif isinstance(edge_attr, str):
                cols = edges_df.columns.intersection([source_col_name, target_col_name, edge_attr])
                edges_df = edges_df[cols]
                use_attrs = edge_attr
            else:
                use_attrs = False
            G = nx.from_pandas_edgelist(edges_df, source=source_col_name, target=target_col_name,
                                        edge_attr=use_attrs, create_using=nx.DiGraph() if directed else nx.Graph())

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
                         target_col_name=target_col_name, edge_attr=edge_attr, filters=filters, directed=directed,
                         relabel_nodes=relabel_nodes, **kwargs)

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
        super().__init__(path, file_resources, source_col_name, target_col_name, edge_attr, filters, directed,
                         relabel_nodes, blocksize, **kwargs)


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
                         target_col_name=target_col_name, edge_attr=edge_attr, filters=filters, directed=directed,
                         relabel_nodes=relabel_nodes, **kwargs)

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        args = dict(na_values=["-"], header=0, low_memory=True,
                    # usecols=['Official Symbol Interactor A', 'Official Symbol Interactor B',
                    #          'Organism Interactor A', 'Score', 'Throughput', 'Qualifications',
                    #          'Modification', 'Phenotypes', 'Source Database'],
                    dtype={'Score': 'float', 'Entrez Gene Interactor A': 'category',
                           'Entrez Gene Interactor B': 'category',
                           'BioGRID ID Interactor A': 'category', 'BioGRID ID Interactor B': 'category',
                           'Systematic Name Interactor A': 'category', 'Systematic Name Interactor B': 'category',
                           'Official Symbol Interactor A': 'category', 'Official Symbol Interactor B': 'category',
                           'Pubmed ID': 'str', 'Throughput': 'category', 'Experimental System Type': 'category',
                           'Experimental System': 'category', 'Modification': 'category', 'Source Database': 'category',
                           'Organism Interactor A': 'category', 'Organism Interactor B': 'category'})

        if blocksize:
            edges = dd.read_table(file_resources["BIOGRID-ALL-LATEST.tab2"], blocksize=blocksize, **args, )
        else:
            edges = pd.read_table(file_resources["BIOGRID-ALL-LATEST.tab2"], **args, )

        self.edges = edges

        return edges

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        df = self.edges
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
                         target_col_name=target_col_name, edge_attr=edge_attr, filters=filters, directed=directed,
                         relabel_nodes=relabel_nodes, **kwargs)

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        df = pd.read_excel(self.file_resources["miRTarBase_MTI.xlsx"])
        self.edges = df
        return df

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        df = self.data
        df = filter_rows(df, filters)

        df['miRNA'] = df['miRNA'].str.rstrip('*')

        if self.strip_mirna_name:
            df['miRNA'] = df['miRNA'].str.lower().str.replace("-3p.*|-5p.*", "", regex=True)

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

    def __init__(self, path='https://dianalab.e-ce.uth.gr/downloads/', file_resources=None, strip_mirna_name=False,
                 source_col_name="mirna", target_col_name="geneId",
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
            file_resources["LncBasev2_download.csv"] = os.path.join(path, "lncbase_v2_exp_data.tar.gz")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, edge_attr=edge_attr, filters=filters, directed=directed,
                         relabel_nodes=relabel_nodes)

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


class TarBase(Interactions):
    """

    """

    def __init__(self, path='https://dianalab.e-ce.uth.gr/downloads', file_resources: Dict = None,
                 source_col_name: str = 'mirna', target_col_name: str = 'geneName',
                 edge_attr: List[str] = None, filters: Union[str, Dict[str, Union[str, List[str]]]] = None,
                 directed: bool = True, relabel_nodes: dict = None, blocksize=None, **kwargs):
        """

        Args:
            path ():
            file_resources ():
            source_col_name ():
            target_col_name ():
            edge_attr ():
            filters ():
            directed ():
            relabel_nodes ():
            blocksize ():
            **kwargs ():
        """
        if file_resources is None:
            file_resources = {
                'tarbase_v8_data.tar.gz': 'https://dianalab.e-ce.uth.gr/downloads/tarbase_v8_data.tar.gz',
                'speclist': 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist',
            }

        super().__init__(path, file_resources, source_col_name, target_col_name, edge_attr, filters, directed,
                         relabel_nodes, blocksize, **kwargs)

    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        edges = pd.read_table(file_resources['tarbase_v8_data.tar.gz'], compression='tar',
                              dtype={'tissue': 'category', 'method': 'category', 'positive_negative': 'category',
                                     'species': 'category',
                                     'direct_indirect': 'category', 'up_down': 'category', 'cell_line': 'category',
                                     })

        if 'speclist' in file_resources:
            species_df = UniProt.get_species_list(file_resources['speclist'])
            species_df = species_df[['Official (scientific) name', 'Common name', 'Synonym']].melt(ignore_index=False)
            species_df = species_df.dropna().reset_index()
            species_name2id = species_df.set_index('value')['NCBI-taxon'].to_dict()
            edges['species_id'] = edges['species'].map(species_name2id)

        self.edges = edges
        return edges

    def load_network(self, file_resources: Dict, source_col_name: str, target_col_name: str, edge_attr: List[str],
                     directed: bool, filters: Dict[str, Any], blocksize=None):
        df = self.data
        df = filter_rows(df, filters)

        # Remove parenthesis containing 3 letter species name
        df['geneName'] = df['geneName'].str.replace(r'(\(\w{3}\)){1}$', '', regex=True)
        idx = df['geneName'].str.contains('\(')
        df.loc[idx, 'geneName'] = df.loc[idx, 'geneName'].str.replace(r'(\(\d of \d\))', '', regex=True).str.strip()

        idx = df['geneName'].str.contains("\(\w*\)", regex=True)
        df.loc[idx, 'geneName'] = df.loc[idx, 'geneName'].str.extract(r'\((\w*)\)(\w*)')[0]

        idx = df['geneName'].str.contains('\(')
        df.loc[idx, 'geneName'] = df.loc[idx, 'geneName'].str.split('(', expand=True)[0]

        g = nx.from_pandas_edgelist(df, source=source_col_name, target=target_col_name,
                                    edge_attr=edge_attr,
                                    create_using=nx.DiGraph() if directed else nx.Graph())
        return g


class RNAInter(Interactions):
    """

    """

    def __init__(self, path='http://www.rnainter.org/raidMedia/download/', file_resources: Dict = None,
                 source_col_name: str = 'Interactor1.Symbol', target_col_name: str = 'Interactor2.Symbol',
                 edge_attr: List[str] = 'score', filters: Union[str, Dict[str, Union[str, List[str]]]] = None,
                 directed: bool = True, relabel_nodes: dict = None, blocksize=None, **kwargs):
        """

        Args:
            path ():
            file_resources ():
            source_col_name ():
            target_col_name ():
            edge_attr ():
            filters ():
            directed ():
            relabel_nodes ():
            blocksize ():
            **kwargs ():
        """
        if file_resources is None:
            file_resources = {
                'Download_data_RR.tar.gz': 'Download_data_RR.tar.gz',
                'Download_data_RP.tar.gz': 'Download_data_RP.tar.gz',
            }

        super().__init__(path, file_resources, source_col_name, target_col_name, edge_attr, filters, directed,
                         relabel_nodes, blocksize, **kwargs)

    def load_dataframe(self, file_resources: Dict, blocksize: int = None) -> pd.DataFrame:
        args = dict(dtype={'Category1': 'category', 'Category2': 'category',
                           'Species1': 'category', 'Species2': 'category', 'score': 'float',
                           'predict': 'category', 'weak': 'category', 'strong': 'category'})
        edge_files = (fn for fn in file_resources if fn.startswith('Download_data'))
        for fn in edge_files:
            if blocksize:
                if not isinstance(file_resources[fn], str): continue
                edges = dd.read_table(file_resources[fn], compression='tar' if fn.endswith('.tar.gz') else None, **args)
            else:
                edges = pd.read_table(file_resources[fn], compression='tar' if fn.endswith('.tar.gz') else None, **args)

        edges = filter_rows(edges, self.filters)

        self.edges = edges
        return edges

    def load_network(self, file_resources, source_col_name='Interactor1.Symbol', target_col_name='Interactor2.Symbol',
                     edge_attr='score', directed=True, filters=None, blocksize=None):
        edges = self.data
        if filters != self.filters:
            edges = filter_rows(edges, filters)

        g = nx.from_pandas_edgelist(edges, source=source_col_name, target=target_col_name,
                                    edge_attr=edge_attr,
                                    create_using=nx.DiGraph() if directed else nx.Graph())
        return g


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
                 edge_attr=["tissue", "positive_negative"], directed=True, relabel_nodes=None, species_id=None,
                 strip_mirna_name=False, **kwargs):
        self.strip_mirna_name = strip_mirna_name
        self.species_id = species_id
        if file_resources is None:
            file_resources = {}
            file_resources["miR_Family_Info.txt.zip"] = os.path.join(path, "miR_Family_Info.txt.zip")
            file_resources["Predicted_Targets_Info.default_predictions.txt"] = os.path.join(path,
                                                                                            "Predicted_Targets_Info.default_predictions.txt")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name,
                         directed=directed, relabel_nodes=relabel_nodes, edge_attr=edge_attr, **kwargs)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        self.df = self.process_miR_family_info_table(file_resources, self.species_id)
        interactions_df = self.process_interactions_table(file_resources, self.df, self.species_id)
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

    def process_interactions_table(self, file_resources, family_to_miR_df, species_id):
        """
        This functions joins the interactions data table between miR Family and targets, and
        Args:
            file_resources:
            family_to_miR_df:
            species_id:

        Returns:

        """
        # Load data frame from file
        family_interactions_df = pd.read_table(file_resources["Predicted_Targets_Info.default_predictions.txt"],
                                               dtype={'Species ID': 'category'},
                                               delimiter='\t', low_memory=True)

        # Select only miRNA-target pairs of certain species_id
        if species_id:
            family_interactions_df = family_interactions_df[family_interactions_df["Species ID"] == species_id]

        family_interactions_df = family_interactions_df.filter(items=["miR Family", "Gene Symbol"], axis="columns")
        family_to_miR_df = family_to_miR_df.filter(items=['miR family', 'MiRBase ID'], axis="columns")
        family_to_miR_df = family_to_miR_df.rename(columns={'miR family': 'miR Family'})

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
                         target_col_name=target_col_name, edge_attr=edge_attr, filters=filters, directed=directed,
                         relabel_nodes=relabel_nodes, verbose=verbose)

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
                 edge_attr=None, filters=None,
                 directed=True, relabel_nodes=None, **kwargs):
        if edge_attr is None:
            edge_attr = ["Interaction Class", "Interaction Mode", "Tissue", "Phenotype"]
        if file_resources is None:
            file_resources = {}
            file_resources["human_interactions.txt"] = os.path.join(path, "human_interactions.txt")

        super().__init__(path, file_resources, source_col_name=source_col_name, target_col_name=target_col_name,
                         edge_attr=edge_attr, filters=filters, directed=directed, relabel_nodes=relabel_nodes, **kwargs)

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

    def __init__(self, path="http://123.59.132.21/lncrna2target/data/", file_resources=None, edge_attr=None,
                 filters=None,
                 directed=True, relabel_nodes=None, version="high_throughput", **kwargs):
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
            file_resources["lncRNA_target_from_high_throughput_experiments.txt.rar"] = \
                os.path.join(path, "lncrna_target.rar")
            file_resources["lncRNA_target_from_low_throughput_experiments.xlsx"] = \
                os.path.join(path, "lncRNA_target_from_low_throughput_experiments.xlsx")

        if self.version == "high_throughput":
            super().__init__(path, file_resources, source_col_name="lncrna_symbol", target_col_name="gene_symbol",
                             edge_attr=edge_attr, filters=filters, directed=directed, relabel_nodes=relabel_nodes,
                             **kwargs)
        if self.version == "low_throughput":
            super().__init__(path, file_resources, source_col_name="GENCODE_gene_name",
                             target_col_name="Target_official_symbol", edge_attr=edge_attr, filters=filters,
                             directed=directed, relabel_nodes=relabel_nodes, **kwargs)

    def load_network(self, file_resources, source_col_name, target_col_name, edge_attr, directed, filters,
                     blocksize=None):
        network = None
        if self.version == "high_throughput":
            network = self.load_network_high_throughput(file_resources, source_col_name, target_col_name, edge_attr,
                                                        directed)
        elif self.version == "low_throughput":
            network = self.load_network_low_throughput(file_resources, source_col_name, target_col_name, edge_attr,
                                                       directed)
        else:
            logger.warn("LncRNA2Target version argument must be one of 'high_throughput' or 'low_throughput'")

        return network

    def load_network_high_throughput(self, file_resources, source_col_name="lncrna_symbol",
                                     target_col_name="gene_symbol",
                                     edge_attr=None, directed=True, filters=None):
        edges = pd.read_table(file_resources["lncRNA_target_from_high_throughput_experiments.txt"], sep="\t")
        edges = filter_rows(edges, filters)

        edges["lncrna_symbol"] = edges["lncrna_symbol"].str.upper()
        edges["lncrna_symbol"] = edges["lncrna_symbol"].str.replace("LINC", "")
        edges["gene_symbol"] = edges["gene_symbol"].str.upper()

        self.data = self.edges = edges
        lncrna2target_high_throughput_network = nx.from_pandas_edgelist(edges,
                                                                        source=source_col_name,
                                                                        target=target_col_name,
                                                                        edge_attr=edge_attr,
                                                                        create_using=nx.DiGraph() if directed else nx.Graph())
        return lncrna2target_high_throughput_network

    def load_network_low_throughput(self, file_resources, source_col_name="GENCODE_gene_name",
                                    target_col_name="Target_official_symbol",
                                    edge_attr=None, directed=True, filters=None):
        edges = pd.read_excel(file_resources["lncRNA_target_from_low_throughput_experiments.xlsx"])
        edges = filter_rows(edges, filters)

        edges["Target_official_symbol"] = edges["Target_official_symbol"].str.replace("(?i)(mir)", "hsa-mir-",
                                                                                      regex=True)
        edges["Target_official_symbol"] = edges["Target_official_symbol"].str.replace("--", "-")
        edges["Target_official_symbol"].apply(lambda x: x.lower() if "mir" in x.lower() else x.upper())
        edges["GENCODE_gene_name"] = edges["GENCODE_gene_name"].str.upper()

        self.data = self.edges = edges
        lncrna2target_low_throughput_network = nx.from_pandas_edgelist(edges,
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
                 edge_attr=["miRNA Interaction Site", "Transcript ID"], directed=True, relabel_nodes=None,
                 **kwargs):
        if file_resources is None:
            file_resources = {}
            file_resources["miRNA_binding_sites.txt"] = os.path.join(path, "miRNA_binding_sites.txt")
            file_resources["general_information.txt"] = os.path.join(path, "general_information.txt")

        super().__init__(path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name,
                         directed=directed, relabel_nodes=relabel_nodes, edge_attr=edge_attr, **kwargs)

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
                 edge_attr=["tarType", "tissueOrCell", "tag", 'class', "level"],
                 filters=None,
                 directed=True, relabel_nodes=None, verbose=False):
        if file_resources is None:
            file_resources = {}
            file_resources["interaction_NPInterv4.expr.txt.gz"] = \
                os.path.join(path, "file/interaction_NPInterv4.expr.txt.gz")

        super().__init__(path=path, file_resources=file_resources, source_col_name=source_col_name,
                         target_col_name=target_col_name, edge_attr=edge_attr, filters=filters, directed=directed,
                         relabel_nodes=relabel_nodes, verbose=verbose)

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
                 min_interactionNum=1, min_expNum=1,
                 edge_attr=None, directed=True, relabel_nodes=None, **kwargs):
        if file_resources is None:
            file_resources = {}
            file_resources["starbase_3.0_lncrna_rna_interactions.csv"] = \
                os.path.join(path, "starbase_3.0_lncrna_rna_interactions.csv")
        self.min_interactionNum = min_interactionNum
        self.min_expNum = min_expNum
        super().__init__(path, file_resources, source_col_name=source_col_name, target_col_name=target_col_name,
                         directed=directed, relabel_nodes=relabel_nodes, edge_attr=edge_attr, **kwargs)

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

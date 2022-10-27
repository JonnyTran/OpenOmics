import os
import warnings
from collections.abc import Iterable
from io import TextIOWrapper, StringIO
from typing import Tuple, List, Dict, Union, Callable, Optional

import dask.dataframe as dd
import networkx as nx
import numpy as np
import obonet
import pandas as pd
import scipy.sparse as ssp
from logzero import logger
from networkx import NetworkXError
from pandas import DataFrame

from openomics.io.read_gaf import read_gaf
from openomics.transforms.adj import slice_adj
from openomics.transforms.agg import get_agg_func
from .base import Database

__all__ = ['GeneOntology', 'UniProtGOA', 'InterPro', 'HumanPhenotypeOntology', ]


class Ontology(Database):
    annotations: pd.DataFrame

    def __init__(self,
                 path,
                 file_resources=None,
                 **kwargs):
        """
        Manages dataset input processing from tables and construct an ontology network from .obo file. There ontology
        network is G(V,E) where there exists e_ij for child i to parent j to present "node i is_a node j".

        Args:
            path:
            file_resources:
            col_rename:
            blocksize:
            verbose:
        """
        # TODO redo the order of .load_datframe(), .load_annotations() and .load_network()
        self.network, self.node_list = self.load_network(file_resources)
        super().__init__(path=path, file_resources=file_resources, **kwargs)

        self.close()

    def load_network(self, file_resources) -> Tuple[nx.MultiDiGraph, List[str]]:
        raise NotImplementedError()

    def filter_network(self, namespace) -> None:
        """
        Filter the subgraph node_list to only `namespace` terms.
        Args:
            namespace: one of {"biological_process", "cellular_component", "molecular_function"}
        """
        terms = self.data[self.data["namespace"] == namespace]["go_id"].unique()
        print("{} terms: {}".format(namespace,
                                    len(terms))) if self.verbose else None
        self.network = self.network.subgraph(nodes=list(terms))
        self.node_list = np.array(list(terms))

    def adj(self, node_list):
        adj_mtx = nx.adj_matrix(self.network, nodelist=node_list)

        if node_list is None or list(node_list) == list(self.node_list):
            return adj_mtx
        elif set(node_list) < set(self.node_list):
            return slice_adj(adj_mtx, list(self.node_list), node_list,
                             None)
        elif not (set(node_list) < set(self.node_list)):
            raise Exception("A node in node_list is not in self.node_list.")

        return adj_mtx

    def filter_annotation(self, annotation: pd.Series):
        go_terms = set(self.node_list)
        filtered_annotation = annotation.map(lambda x: list(set(x) & go_terms)
                                             if isinstance(x, list) else [])

        return filtered_annotation

    def get_child_nodes(self):
        adj = self.adj(self.node_list)
        leaf_terms = self.node_list[np.nonzero(adj.sum(axis=0) == 0)[1]]
        return leaf_terms

    def get_root_nodes(self):
        adj = self.adj(self.node_list)
        parent_terms = self.node_list[np.nonzero(adj.sum(axis=1) == 0)[0]]
        return parent_terms

    def get_dfs_paths(self, root_nodes: list, filter_duplicates=False):
        """
        Return all depth-first search paths from root node(s) to children node by traversing the ontology directed graph.
        Args:
            root_nodes (list): ["GO:0008150"] if biological processes, ["GO:0003674"] if molecular_function, or ["GO:0005575"] if cellular_component
            filter_duplicates (bool): whether to remove duplicated paths that end up at the same leaf nodes

        Returns: pd.DataFrame of all paths starting from the root nodes.
        """
        if not isinstance(root_nodes, list):
            root_nodes = list(root_nodes)

        paths = list(dfs_path(self.network, root_nodes))
        paths = list(flatten_list(paths))
        paths_df = pd.DataFrame(paths)

        if filter_duplicates:
            paths_df = paths_df[~paths_df.duplicated(keep="first")]
            paths_df = filter_dfs_paths(paths_df)

        return paths_df

    def remove_predecessor_terms(self, annotation: pd.Series, sep="\||;"):
        # leaf_terms = self.get_child_nodes()
        # if not annotation.map(lambda x: isinstance(x, (list, np.ndarray))).any() and sep:
        #     annotation = annotation.str.split(sep)
        #
        # parent_terms = annotation.map(lambda x: list(
        #     set(x) & set(leaf_terms)) if isinstance(x, (list, np.ndarray)) else None)
        # return parent_terms
        raise NotImplementedError

    def get_subgraph(self, edge_types: Union[str, List[str]]) -> Union[nx.MultiDiGraph, nx.DiGraph]:
        if not hasattr(self, "_subgraphs"):
            self._subgraphs = {}
        elif edge_types in self._subgraphs:
            return self._subgraphs[edge_types]

        if edge_types and isinstance(self.network, (nx.MultiGraph, nx.MultiDiGraph)):
            # Needed to create new nx.Graph because .edge_subgraph is too slow to iterate on (idk why)
            g = nx.from_edgelist([(u, v) for u, v, k in self.network.edges if k in edge_types],
                                 create_using=nx.DiGraph if self.network.is_directed() else nx.Graph)
        else:
            raise Exception("Must provide `edge_types` keys for a nx.MultiGraph type.")

        self._subgraphs[edge_types] = g

        return g

    def add_predecessor_terms(self, anns: pd.Series, edge_type: Union[str, List[str]] = 'is_a', sep="\||;"):
        anns_w_parents = anns.map(lambda x: [] if not isinstance(x, (list, np.ndarray)) else x) + \
                         get_predecessor_terms(anns, self.get_subgraph(edge_type))

        return anns_w_parents

    @staticmethod
    def get_node_color(file="~/Bioinformatics_ExternalData/GeneOntology/go_colors_biological.csv", ):
        go_colors = pd.read_csv(file)

        def selectgo(x):
            terms = [term for term in x if isinstance(term, str)]
            if len(terms) > 0:
                return terms[-1]
            else:
                return None

        go_colors["node"] = go_colors[[
            col for col in go_colors.columns if col.isdigit()
        ]].apply(selectgo, axis=1)
        go_id_colors = go_colors[go_colors["node"].notnull()].set_index("node")["HCL.color"]
        go_id_colors = go_id_colors[~go_id_colors.index.duplicated(keep="first")]

        print(go_id_colors.unique().shape, go_colors["HCL.color"].unique().shape)
        return go_id_colors

    def split_annotations(self, src_node_col="gene_name", dst_node_col="go_id", groupby: List[str] = ["Qualifier"],
                          train_date="2017-06-15", valid_date="2017-11-15", test_date="2021-12-31",
                          query: Optional[str] = "Evidence in ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC']",
                          filter_src_nodes: pd.Index = None, filter_dst_nodes: pd.Index = None,
                          agg: Union[Callable, str] = "unique") -> Tuple[DataFrame, DataFrame, DataFrame]:
        """

        Args:
            src_node_col (str): Name of column containg the the src node types.
            dst_node_col (str): Name of column containg the the dst node types.
            train_date (str): A date before which the annotations belongs in the training set.
            valid_date (str): A date before which the annotations belongs in the validation set.
            test_date (str): A date before which the annotations belongs in the testing set.
            groupby (str): A list of strings to groupby annotations on, default [`src_node_col`, "Qualifier"].
            query (str, optional): A pandas query string to filter annotations. Default, only select ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC'] annotations.
            filter_src_nodes (pd.Index): A subset annotations by these values on `src_node_col`.
            filter_dst_nodes (pd.Index): A subset annotations by these values on `dst_node_col`.
            agg (str): Either "unique" or "add_parent", or a callable function, or a dd.Aggregation() for aggregating on the `dst_node_col` column after groupby on `groupby`.
        """
        raise NotImplementedError


class GeneOntology(Ontology):
    """Loads the GeneOntology database from http://geneontology.org .

    Default path: "http://geneontology.org/gene-associations/".

    Default file_resources: {
        "go-basic.obo": "http://purl.obolibrary.org/obo/go/go-basic.obo",
        "goa_human.gaf": "goa_human.gaf.gz",
        "goa_human_rna.gaf": "goa_human_rna.gaf.gz",
        "goa_human_isoform.gaf": "goa_human_isoform.gaf.gz",
    }
    """
    COLUMNS_RENAME_DICT = {
        "DB_Object_Symbol": "gene_name",
        "DB_Object_ID": "gene_id",
        "GO_ID": "go_id",
        "Taxon_ID": 'species_id',
    }

    DROP_COLS = {'DB:Reference', 'With', 'Annotation_Extension', 'Gene_Product_Form_ID'}

    def __init__(
        self,
        path="http://geneontology.org/gene-associations/",
        species='human',
        file_resources=None,
        index_col='DB_Object_Symbol',
        keys=None,
        col_rename=COLUMNS_RENAME_DICT,
        blocksize=None,
        **kwargs
    ):
        """
        Loads the GeneOntology database from http://geneontology.org .

            Default path: "http://geneontology.org/gene-associations/" .
            Default file_resources: {
                "go-basic.obo": "http://purl.obolibrary.org/obo/go/go-basic.obo",
                "goa_human.gaf": "goa_human.gaf.gz",
                "goa_human_rna.gaf": "goa_human_rna.gaf.gz",
                "goa_human_isoform.gaf": "goa_human_isoform.gaf.gz",
            }

        Data for GO term annotations in .gpi files are already included in .obo file, so this module doesn't maker use of .gpi files.

        Handles downloading the latest Gene Ontology obo and annotation data, preprocesses them. It provide
        functionalities to create a directed acyclic graph of GO terms, filter terms, and filter annotations.
        """
        if species and not hasattr(self, 'species'):
            self.species = species.lower()
        elif species is None:
            self.species = 'uniprot'

        if file_resources is None:
            file_resources = {
                f"goa_{self.species}.gaf.gz": f"goa_{self.species}.gaf.gz",
            }
            if species != 'uniprot':
                file_resources[f"goa_{self.species}_rna.gaf.gz"] = f"goa_{self.species}_rna.gaf.gz"
                file_resources[f"goa_{self.species}_isoform.gaf.gz"] = f"goa_{self.species}_isoform.gaf.gz"

        if not any('.obo' in file for file in file_resources):
            warnings.warn(
                f'No .obo file provided in `file_resources`, so automatically adding "http://purl.obolibrary.org/obo/go/go-basic.obo"')
            file_resources["go-basic.obo"] = "http://purl.obolibrary.org/obo/go/go-basic.obo"

        super().__init__(path, file_resources, index_col=index_col, keys=keys, col_rename=col_rename,
                         blocksize=blocksize,
                         **kwargs)

    def info(self):
        print("network {}".format(nx.info(self.network)))

    def load_dataframe(self, file_resources: Dict[str, TextIOWrapper], blocksize=None) -> DataFrame:
        if self.network:
            # Annotations for each GO term from nodes in the NetworkX graph created by the .obo file
            go_terms = pd.DataFrame.from_dict(dict(self.network.nodes(data=True)), orient='index')
            go_terms["def"] = go_terms["def"].apply(
                lambda x: x.split('"')[1] if isinstance(x, str) else None)
            go_terms.index.name = "go_id"
        else:
            go_terms = None

        # Handle .gaf annotation files
        dfs = {}
        for filename, filepath_or_buffer in file_resources.items():
            gaf_name = filename.split(".")[0]
            # Ensure no duplicate GAF file (if having files uncompressed with same prefix)
            if gaf_name in dfs: continue

            if blocksize and isinstance(filepath_or_buffer, str):
                if filename.endswith(".processed.parquet"):
                    # Parsed and filtered gaf file
                    dfs[gaf_name] = dd.read_parquet(filepath_or_buffer, chunksize=blocksize)
                    if dfs[gaf_name].index.name != self.index_col and self.index_col in dfs[gaf_name].columns:
                        dfs[gaf_name] = dfs[gaf_name].set_index(self.index_col, sorted=True)
                    if not dfs[gaf_name].known_divisions:
                        dfs[gaf_name].divisions = dfs[gaf_name].compute_current_divisions()

                elif (filename.endswith(".parquet") or filename.endswith(".gaf")):
                    # .parquet from .gaf.gz file, unfiltered, with raw str values
                    dfs[gaf_name] = read_gaf(filepath_or_buffer, blocksize=blocksize, index_col=self.index_col,
                                             keys=self.keys, usecols=self.usecols)

                elif filename.endswith(".gaf.gz"):
                    # Compressed .gaf file downloaded
                    dfs[gaf_name] = read_gaf(filepath_or_buffer, blocksize=blocksize, index_col=self.index_col,
                                             keys=self.keys, usecols=self.usecols, compression='gzip')

            else:
                if filename.endswith(".processed.parquet"):
                    dfs[gaf_name] = pd.read_parquet(filepath_or_buffer)
                if filename.endswith(".gaf"):
                    dfs[gaf_name] = read_gaf(filepath_or_buffer, index_col=self.index_col, keys=self.keys,
                                             usecols=self.usecols)

        # Filter and set index divisions
        # for gaf_name, df in dfs.items():
        #     if self.keys is not None and self.index_col and df.index.name == self.index_col:
        #         dfs[gaf_name] = df.loc[df.index.isin(self.keys)]
        #     elif self.keys is not None and self.index_col and df.index.name != self.index_col:
        #         dfs[gaf_name] = df.loc[df[self.index_col].isin(self.keys)]
        #
        #     if isinstance(df, dd.DataFrame) and not df.known_divisions:
        #         dfs[gaf_name].divisions = df.compute_current_divisions()

        if len(dfs):
            self.annotations = dd.concat(list(dfs.values()), interleave_partitions=True) \
                if blocksize else pd.concat(dfs.values())

            if len(self.annotations.columns.intersection(UniProtGOA.COLUMNS_RENAME_DICT.keys())):
                self.annotations = self.annotations.rename(columns=UniProtGOA.COLUMNS_RENAME_DICT)
                if self.annotations.index.name in UniProtGOA.COLUMNS_RENAME_DICT:
                    self.annotations.index = self.annotations.index.rename(
                        UniProtGOA.COLUMNS_RENAME_DICT[self.annotations.index.name])

        return go_terms

    def load_network(self, file_resources) -> Tuple[nx.Graph, np.ndarray]:
        for file in file_resources:
            if file.endswith(".obo"):
                network: nx.MultiDiGraph = obonet.read_obo(file_resources[file])
                network = network.reverse(copy=True)
                node_list = np.array(network.nodes)

                return network, node_list

    def split_annotations(self, src_node_col="gene_name", dst_node_col="go_id", groupby: List[str] = ["Qualifier"],
                          train_date="2017-06-15", valid_date="2017-11-15", test_date="2021-12-31",
                          query: str = "Evidence in ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC']",
                          filter_src_nodes: pd.Index = None, filter_dst_nodes: pd.Index = None,
                          agg: Union[str, Callable, dd.Aggregation] = "unique") \
        -> Tuple[DataFrame, DataFrame, DataFrame]:
        assert isinstance(groupby, list) and groupby, f"`groupby` must be a nonempty list of strings. Got {groupby}"

        # Set the source column (i.e. protein_id or gene_name), to be the first in groupby
        if src_node_col not in groupby:
            groupby = [src_node_col] + groupby
        if "Qualifier" not in groupby and "Qualifier" in self.annotations.columns:
            groupby.append("Qualifier")

        # Aggregator function
        if agg == "add_parent":
            subgraph = self.get_subgraph(edge_types="is_a")
            node_ancestors = {node: nx.ancestors(subgraph, node) for node in subgraph.nodes}

            if isinstance(self.annotations, dd.DataFrame):
                agg = dd.Aggregation(name='_unique_add_parent',
                                     chunk=lambda s: s.unique(),
                                     agg=lambda s0: s0.apply(get_predecessor_terms, node_ancestors, keep_terms=True),
                                     finalize=lambda s1: s1.apply(lambda li: np.hstack(li) if li else None))
            else:
                agg = lambda s: get_predecessor_terms(s, g=node_ancestors, join_groups=True, keep_terms=True)

        elif agg == 'unique' and isinstance(self.annotations, dd.DataFrame):
            agg = get_agg_func('unique', use_dask=True)

        elif isinstance(self.annotations, dd.DataFrame) and not isinstance(agg, dd.Aggregation):
            raise Exception("`agg` must be a dd.Aggregation for groupby.agg() on columns of a dask DataFrame")

        def _remove_dup_neg_go_id(s: pd.Series) -> pd.Series:
            if s.isna().any():
                return s
            elif isinstance(s[neg_dst_col], Iterable) and isinstance(s[dst_node_col], Iterable):
                rm_dups_go_id = [go_id for go_id in s[neg_dst_col] if go_id not in s[dst_node_col]]
                if len(rm_dups_go_id) == 0:
                    rm_dups_go_id = None
                s[neg_dst_col] = rm_dups_go_id
            return s

        neg_dst_col = f"neg_{dst_node_col}"

        # Filter annotations
        annotations = self.annotations
        if query:
            annotations = annotations.query(query)
        if filter_src_nodes is not None:
            annotations = annotations[annotations[src_node_col].isin(filter_src_nodes)]
        if filter_dst_nodes is not None:
            annotations = annotations[annotations[dst_node_col].isin(filter_dst_nodes)]
        if annotations.index.name in groupby:
            annotations = annotations.reset_index()

        # Split train/valid/test annotations
        train_anns = annotations[annotations["Date"] <= pd.to_datetime(train_date)]
        valid_anns = annotations[(annotations["Date"] <= pd.to_datetime(valid_date)) & \
                                 (annotations["Date"] > pd.to_datetime(train_date))]
        test_anns = annotations[(annotations["Date"] <= pd.to_datetime(test_date)) & \
                                (annotations["Date"] > pd.to_datetime(valid_date))]

        outputs = []
        for anns in [train_anns, valid_anns, test_anns]:
            # Keep track of which annotation has a "NOT" Qualifier
            is_neg_ann = anns["Qualifier"].map(lambda li: "NOT" in li)

            # Convert `Qualifiers` entries of list of strings to string
            args = dict(meta=pd.Series([""])) if isinstance(anns, dd.DataFrame) else {}
            anns.loc[:, 'Qualifier'] = anns['Qualifier'].apply(
                lambda li: "".join([i for i in li if i != "NOT"]), **args)

            # Aggregate gene-GO annotations
            if isinstance(anns, pd.DataFrame) and len(anns.index):
                pos_anns = anns[~is_neg_ann].groupby(groupby).agg({dst_node_col: agg})
                neg_anns = anns[is_neg_ann].groupby(groupby).agg(**{neg_dst_col: (dst_node_col, agg)})
                pos_neg_anns = pd.concat([pos_anns, neg_anns], axis=1)

            elif isinstance(anns, dd.DataFrame) and len(anns.index) and dst_node_col in anns.columns:
                pos_anns = anns[~is_neg_ann].groupby(groupby).agg({dst_node_col: agg})
                if False and len(is_neg_ann.index):
                    neg_anns = anns[is_neg_ann].groupby(groupby).agg({dst_node_col: agg})
                    neg_anns.columns = [neg_dst_col]
                    pos_neg_anns = dd.concat([pos_anns, neg_anns], axis=1)
                else:
                    pos_neg_anns = pos_anns
                    pos_neg_anns[neg_dst_col] = None

            else:
                pos_neg_anns = pd.DataFrame(
                    columns=[dst_node_col, neg_dst_col],
                    index=pd.MultiIndex(levels=[[] for i in range(len(groupby))],
                                        codes=[[] for i in range(len(groupby))], names=groupby))
                outputs.append(pos_neg_anns)
                continue

            if isinstance(pos_neg_anns, pd.DataFrame):
                pos_neg_anns = pos_neg_anns.drop([""], axis='index', errors="ignore")

            # Remove "GO:0005515" (protein binding) annotations for a gene if it's the gene's only annotation
            _exclude_single_fn = lambda li: None \
                if isinstance(li, Iterable) and len(li) == 1 and "GO:0005515" in li else li
            args = dict(meta=pd.Series([list()])) if isinstance(anns, dd.DataFrame) else {}
            pos_neg_anns.loc[:, dst_node_col] = pos_neg_anns[dst_node_col].apply(_exclude_single_fn, **args)

            # Drop rows with all nan values
            if isinstance(pos_neg_anns, pd.DataFrame):
                pos_neg_anns = pos_neg_anns.drop(pos_neg_anns.index[pos_neg_anns.isna().all(1)], axis='index')

            # Ensure no negative terms duplicates positive annotations
            if len(is_neg_ann.index):
                args = dict(meta=pd.DataFrame({dst_node_col: [], neg_dst_col: []})) \
                    if isinstance(anns, dd.DataFrame) else {}
                pos_neg_anns = pos_neg_anns.apply(_remove_dup_neg_go_id, axis=1, **args)

            outputs.append(pos_neg_anns)

        return tuple(outputs)


def get_predecessor_terms(anns: Union[pd.Series, Iterable], g: Union[Dict[str, List[str]], nx.MultiDiGraph],
                          join_groups=False, keep_terms=True, exclude={'GO:0005575', 'GO:0008150', 'GO:0003674'}) \
    -> Union[pd.Series, List[str]]:
    """

    Args:
        anns ():
        g (nx.MultiDiGraph, Dict[str,Set[str]]): Either a NetworkX DAG or a precomputed lookup table of node to ancestors
        join_groups (): whether to concatenate multiple
        keep_terms ():
        exclude ():

    Returns:

    """
    if exclude is None:
        exclude = {}

    def _get_ancestors(terms: Iterable):
        try:
            if isinstance(terms, Iterable):
                if isinstance(g, dict):
                    parents = {parent \
                               for term in terms if term in g \
                               for parent in g[term] if parent not in exclude}

                elif isinstance(g, nx.Graph):
                    parents = {parent \
                               for term in terms if term in g.nodes \
                               for parent in nx.ancestors(g, term) if parent not in exclude}
                else:
                    raise Exception("Provided `g` arg must be either an nx.Graph or a Dict")
            else:
                parents = []

            if keep_terms and isinstance(terms, list):
                terms.extend(parents)
                out = terms
            else:
                out = list(parents)

        except (NetworkXError, KeyError) as nxe:
            if "foo" not in nxe.__str__():
                logger.error(f"{nxe.__class__.__name__} get_predecessor_terms._get_ancestors: {nxe}")
            out = terms if keep_terms else []

        return out

    if isinstance(anns, pd.Series):
        if (anns.map(type) == str).any():
            anns = anns.map(lambda s: [s])

        parent_terms = anns.map(_get_ancestors)
        if join_groups:
            parent_terms = sum(parent_terms, [])

    elif isinstance(anns, Iterable):
        parent_terms = _get_ancestors(anns)

    else:
        parent_terms = []

    return parent_terms

class UniProtGOA(GeneOntology):
    """Loads the GeneOntology database from https://www.ebi.ac.uk/GOA/ .

    Default path: "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/" .
    Default file_resources: {
        "goa_uniprot_all.gaf": "goa_uniprot_all.gaf.gz",
    }
    """

    COLUMNS_RENAME_DICT = {
        "DB_Object_ID": "protein_id",
        "DB_Object_Symbol": "gene_name",
        "GO_ID": "go_id",
        "Taxon_ID": 'species_id',
    }
    def __init__(
        self,
        path="ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/",
        species="HUMAN",
        file_resources=None,
        index_col='DB_Object_ID', keys=None,
        col_rename=COLUMNS_RENAME_DICT,
        blocksize=None,
        **kwargs,
    ):
        """
        Loads the UniProtGOA database from https://www.ebi.ac.uk/GOA/ .

            Default path: "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/" .
            Default file_resources: {
                "goa_uniprot_all.gaf.gz": "goa_uniprot_all.gaf.gz",
                "go.obo": "http://current.geneontology.org/ontology/go.obo",
            }

        Handles downloading the latest Gene Ontology obo and annotation data, preprocesses them. It provides
        functionalities to create a directed acyclic graph of GO terms, filter terms, and filter annotations.

        Args:
            path ():
            species ():
            file_resources ():
            index_col ():
            keys ():
            col_rename ():
            blocksize ():
            **kwargs ():
        """
        if species is None:
            self.species = species = 'UNIPROT'
            substr = 'uniprot_all'
        else:
            self.species = species.upper()
            substr = species.lower()

        if file_resources is None:
            file_resources = {
                "go.obo": "http://current.geneontology.org/ontology/go.obo",
                f"goa_{self.species.lower()}.gaf.gz": os.path.join(species, f"goa_{substr}.gaf.gz"),
                # f"goa_{self.species.lower()}_isoform.gaf.gz": os.path.join(species, f"goa_{substr}_isoform.gaf.gz"),
                # f"goa_{self.species.lower()}_complex.gaf.gz": os.path.join(species, f"goa_{substr}_complex.gaf.gz"),
            }

        if not any('.obo' in file for file in file_resources):
            warnings.warn(f'No .obo file provided in `file_resources`, '
                          f'so automatically adding "http://purl.obolibrary.org/obo/go/go-basic.obo"')
            file_resources["go-basic.obo"] = "http://purl.obolibrary.org/obo/go/go-basic.obo"

        super().__init__(path=path, file_resources=file_resources, index_col=index_col, keys=keys,
                         col_rename=col_rename,
                         blocksize=blocksize, **kwargs)


class InterPro(Ontology):
    """
    Default parameters
    path="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/"
    file_resources = {}
    file_resources["entry.list"] = os.path.join(path, "entry.list")
    file_resources["protein2ipr.dat.gz"] = os.path.join(path, "protein2ipr.dat.gz")
    file_resources["interpro2go"] = os.path.join(path, "interpro2go")
    file_resources["ParentChildTreeFile.txt"] = os.path.join(path, "ParentChildTreeFile.txt")
    """

    def __init__(self, path="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/", index_col='UniProtKB-AC',
                 keys=None,
                 file_resources=None, col_rename=None, **kwargs):
        """
        Default parameters
            path="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/"
            file_resources = {}
            file_resources["entry.list"] = os.path.join(path, "entry.list")
            file_resources["protein2ipr.dat.gz"] = os.path.join(path, "protein2ipr.dat.gz")
            file_resources["interpro2go"] = os.path.join(path, "interpro2go")
            file_resources["ParentChildTreeFile.txt"] = os.path.join(path, "ParentChildTreeFile.txt")

        Args:
            path:
            file_resources:
            col_rename:
            verbose:
        """
        assert keys is not None
        assert index_col is not None

        if file_resources is None:
            file_resources = {}
            file_resources["entry.list"] = os.path.join(path, "entry.list")
            file_resources["protein2ipr.dat.gz"] = os.path.join(path, "protein2ipr.dat.gz")
            file_resources["interpro2go"] = os.path.join(path, "interpro2go")
            file_resources["ParentChildTreeFile.txt"] = os.path.join(path, "ParentChildTreeFile.txt")

        super().__init__(path=path, file_resources=file_resources, index_col=index_col, keys=keys,
                         col_rename=col_rename, **kwargs)

    def load_dataframe(self, file_resources: Dict[str, TextIOWrapper], blocksize=None):
        ipr_entries = pd.read_table(file_resources["entry.list"], index_col="ENTRY_AC")
        ipr2go = self.parse_interpro2go(file_resources["interpro2go"])
        if ipr2go is not None:
            ipr_entries = ipr_entries.join(ipr2go.groupby('ENTRY_AC')["go_id"].unique(), on="ENTRY_AC")

        # Use Dask
        args = dict(names=['UniProtKB-AC', 'ENTRY_AC', 'ENTRY_NAME', 'accession', 'start', 'stop'],
                    usecols=['UniProtKB-AC', 'ENTRY_AC', 'start', 'stop'],
                    dtype={'UniProtKB-AC': 'category', 'ENTRY_AC': 'category', 'start': 'int8', 'stop': 'int8'},
                    low_memory=True,
                    blocksize=None if isinstance(blocksize, bool) else blocksize)
        if 'protein2ipr.parquet' in file_resources:
            annotations = dd.read_parquet(file_resources["protein2ipr.parquet"])
        else:
            annotations = dd.read_table(file_resources["protein2ipr.dat"], **args)

        if self.keys is not None and self.index_col in annotations.columns:
            annotations = annotations.loc[annotations[self.index_col].isin(self.keys)]
        elif self.keys is not None and self.index_col == annotations.index.name:
            annotations = annotations.loc[annotations.index.isin(self.keys)]

        # if annotations.index.name != self.index_col:
        #     annotations = annotations.set_index(self.index_col, sorted=True)
        # if not annotations.known_divisions:
        #     annotations.divisions = annotations.compute_current_divisions()

        # Set ordering for rows and columns
        row_order = self.keys
        col_order = ipr_entries.index
        row2idx = {node: i for i, node in enumerate(row_order)}
        col2idx = {node: i for i, node in enumerate(col_order)}

        def edgelist2coo(edgelist_df: DataFrame, source='UniProtKB-AC', target='ENTRY_AC') -> Optional[ssp.coo_matrix]:
            if edgelist_df.shape[0] == 1 and edgelist_df.iloc[0, 0] == 'foo':
                return None

            if edgelist_df.index.name == source:
                source_nodes = edgelist_df.index
            else:
                source_nodes = edgelist_df[source]

            edgelist_df = edgelist_df.assign(row=source_nodes.map(row2idx).astype('int'),
                                             col=edgelist_df[target].map(col2idx).astype('int'))

            edgelist_df = edgelist_df.dropna(subset=['row', 'col'])
            if edgelist_df.shape[0] == 0:
                return None

            values = np.ones(edgelist_df.index.size)
            coo = ssp.coo_matrix((values, (edgelist_df['row'], edgelist_df['col'])),
                                 shape=(len(row2idx), ipr_entries.index.size))
            return coo

        # Create a sparse adjacency matrix each partition, then combine them
        adj = annotations.reduction(chunk=edgelist2coo,
                                    aggregate=lambda x: x.dropna().sum() if not x.isna().all() else None,
                                    meta=pd.Series([ssp.coo_matrix])).compute()
        assert len(adj) == 1, f"len(adj) = {len(adj)}"

        # Create a sparse matrix of UniProtKB-AC x ENTRY_AC
        self.annotations = pd.DataFrame.sparse.from_spmatrix(adj[0], index=row_order, columns=col_order)

        return ipr_entries

    def parse_interpro2go(self, file: StringIO) -> pd.DataFrame:
        def _process_line(line: str) -> Tuple[str, str, str]:
            pos = line.find('> GO')
            interpro_terms, go_term = line[:pos], line[pos:]
            interpro_id, interpro_name = interpro_terms.strip().split(' ', 1)
            go_name, go_id = go_term.split(';')
            go_desc = go_name.strip('> GO:')

            return (interpro_id.strip().split(':')[1], go_id.strip(), go_desc)

        if isinstance(file, str):
            with open(os.path.expanduser(file), 'r') as file:
                tuples = [_process_line(line.strip()) for line in file if line[0] != '!']

            ipr2go = pd.DataFrame(tuples, columns=['ENTRY_AC', "go_id", "go_desc"])
            return ipr2go

    def load_network(self, file_resources) -> Tuple[nx.Graph, np.ndarray]:
        network, node_list = None, None
        for filename in file_resources:
            if 'ParentChildTreeFile' in filename and isinstance(file_resources[filename], str):
                network: nx.MultiDiGraph = self.parse_ipr_treefile(file_resources[filename])
                node_list = np.array(network.nodes)

        return network, node_list

    def parse_ipr_treefile(self, lines: Union[List[str], StringIO]) -> nx.MultiDiGraph:
        """Parse the InterPro Tree from the given file.
        Args:
            lines: A readable file or file-like
        """
        if isinstance(lines, str):
            lines = open(os.path.expanduser(lines), 'r')

        graph = nx.MultiDiGraph()
        previous_depth, previous_name = 0, None
        stack = [previous_name]

        def count_front(s: str) -> int:
            """Count the number of leading dashes on a string."""
            for position, element in enumerate(s):
                if element != '-':
                    return position

        for line in lines:
            depth = count_front(line)
            interpro_id, name, *_ = line[depth:].split('::')

            if depth == 0:
                stack.clear()
                stack.append(interpro_id)

                graph.add_node(interpro_id, interpro_id=interpro_id, name=name)

            else:
                if depth > previous_depth:
                    stack.append(previous_name)

                elif depth < previous_depth:
                    del stack[-1]

                parent = stack[-1]

                graph.add_node(interpro_id, interpro_id=interpro_id, parent=parent, name=name)
                graph.add_edge(parent, interpro_id, key="is_a")

            previous_depth, previous_name = depth, interpro_id

        lines.close()
        return graph


class HumanPhenotypeOntology(Ontology):
    """Loads the Human Phenotype Ontology database from https://hpo.jax.org/app/ .

        Default path: "http://geneontology.org/gene-associations/" .
        Default file_resources: {
            "hp.obo": "http://purl.obolibrary.org/obo/hp.obo",
        }
        """

    COLUMNS_RENAME_DICT = {}

    def __init__(
        self,
        path="https://hpo.jax.org/",
        file_resources=None,
        col_rename=COLUMNS_RENAME_DICT,
        blocksize=0,
        verbose=False,
    ):
        """
        Handles downloading the latest Human Phenotype Ontology obo and annotation data, preprocesses them. It provide
        functionalities to create a directed acyclic graph of Ontology terms, filter terms, and filter annotations.
        """
        if file_resources is None:
            file_resources = {
                "hp.obo": "http://purl.obolibrary.org/obo/hp.obo",
            }
        super().__init__(
            path,
            file_resources,
            col_rename=col_rename,
            blocksize=blocksize,
            verbose=verbose,
        )

    def info(self):
        print("network {}".format(nx.info(self.network)))
    def load_network(self, file_resources):
        for file in file_resources:
            if ".obo" in file:
                network = obonet.read_obo(file_resources[file])
                network = network.reverse(copy=True)
                node_list = np.array(network.nodes)
        return network, node_list




def traverse_predecessors(network, seed_node, type=["is_a", "part_of"]):
    """
    Returns all successor terms from seed_node by traversing the ontology network with edges == `type`.
    Args:
        seed_node: seed node of the traversal
        type: the ontology type to include
    Returns:
        generator of list of lists for each dfs branches.
    """
    parents = dict(network.pred[seed_node])
    for parent, v in parents.items():
        if list(v.keys())[0] in type:
            yield [parent] + list(traverse_predecessors(network, parent, type))


def flatten(lst):
    return sum(([x] if not isinstance(x, list) else flatten(x) for x in lst),
               [])


def dfs_path(graph, path):
    node = path[-1]
    successors = list(graph.successors(node))
    if len(successors) > 0:
        for child in successors:
            yield list(dfs_path(graph, path + [child]))
    else:
        yield path


def flatten_list(list_in):
    if isinstance(list_in, list):
        for l in list_in:
            if isinstance(list_in[0], list):
                for y in flatten_list(l):
                    yield y
            elif isinstance(list_in[0], str):
                yield list_in
    else:
        yield list_in


def filter_dfs_paths(paths_df: pd.DataFrame):
    idx = {}
    for col in sorted(paths_df.columns[:-1], reverse=True):
        idx[col] = ~(paths_df[col].notnull()
                     & paths_df[col].duplicated(keep="first")
                     & paths_df[col + 1].isnull())

    idx = pd.DataFrame(idx)

    paths_df = paths_df[idx.all(axis=1)]
    return paths_df


def write_taxonomy(network, root_nodes, file_path):
    """

    Args:
        network: A network with edge(i, j) where i is a node and j is a child of i.
        root_nodes (list): a list of node names
        file_path (str):
    """
    file = open(file_path, "a")
    file.write("Root\t" + "\t".join(root_nodes) + "\n")

    for root_node in root_nodes:
        for node, children in nx.traversal.bfs_successors(network, root_node):
            if len(children) > 0:
                file.write(node + "\t" + "\t".join(children) + "\n")
    file.close()

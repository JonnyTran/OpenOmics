import os
from io import TextIOWrapper, StringIO
from typing import Tuple, List, Dict, Iterable, Union

import dask.dataframe as dd
import networkx as nx
import numpy as np
import obonet
import pandas as pd
import scipy.sparse as ssp
import tqdm
from Bio.UniProt.GOA import _gaf20iterator, _gaf10iterator
from pandas import DataFrame

from .base import Database
from ..utils.df import slice_adj


class Ontology(Database):
    def __init__(self,
                 path,
                 file_resources=None,
                 col_rename=None,
                 npartitions=0,
                 verbose=False):
        """
        Manages dataset input processing from tables and construct an ontology network from .obo file. There ontology
        network is G(V,E) where there exists e_ij for child i to parent j to present "node i is_a node j".

        Args:
            path:
            file_resources:
            col_rename:
            npartitions:
            verbose:
        """
        self.network, self.node_list = self.load_network(file_resources)

        super().__init__(
            path=path,
            file_resources=file_resources,
            col_rename=col_rename,
            npartitions=npartitions,
            verbose=verbose,
        )

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
        leaf_terms = self.get_child_nodes()
        if not annotation.map(lambda x: isinstance(x, list)).any() and sep:
            annotation = annotation.str.split(sep)

        go_terms_parents = annotation.map(lambda x: list(
            set(x) & set(leaf_terms)) if isinstance(x, list) else None)
        return go_terms_parents

    def get_predecessor_terms(self, annotations: pd.Series, edge_type='is_a'):
        if isinstance(self.network, (nx.MultiGraph, nx.MultiDiGraph)):
            ontology_g = self.network.edge_subgraph([(u, v, k) for u, v, k in self.network.edges if k == edge_type])
        else:
            ontology_g = self.network

        go_terms_parents = annotations.map(
            lambda annotations: \
                list({parent for term in annotations \
                      for parent in list(nx.ancestors(ontology_g, term))}) \
                    if isinstance(annotations, list) else [])
        return go_terms_parents

    def add_predecessor_terms(self, annotation: pd.Series, edge_type='is_a', sep="\||;", return_str=False):
        if sep and annotation.dtype == np.object and annotation.str.contains(sep, regex=True).any():
            ann_lists = annotation.str.split(sep)
        else:
            ann_lists = annotation

        ann_with_parents = ann_lists + self.get_predecessor_terms(ann_lists, edge_type)

        if return_str:
            ann_with_parents = ann_with_parents.map(
                lambda x: "|".join(x) if isinstance(x, list) else None)

        return ann_with_parents

    @staticmethod
    def get_node_color(
        file="~/Bioinformatics_ExternalData/GeneOntology/go_colors_biological.csv",
    ):
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
        go_id_colors = go_colors[go_colors["node"].notnull()].set_index(
            "node")["HCL.color"]
        go_id_colors = go_id_colors[~go_id_colors.index.duplicated(
            keep="first")]
        print(go_id_colors.unique().shape,
              go_colors["HCL.color"].unique().shape)
        return go_id_colors



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
        npartitions=0,
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
            npartitions=npartitions,
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


class GeneOntology(Ontology):
    """Loads the GeneOntology database from http://geneontology.org .

    Default path: "http://geneontology.org/gene-associations/" .
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
    }

    def __init__(
        self,
        path="http://geneontology.org/gene-associations/",
        file_resources=None,
        col_rename=COLUMNS_RENAME_DICT,
        npartitions=0,
        verbose=False,
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

        Handles downloading the latest Gene Ontology obo and annotation data, preprocesses them. It provide
        functionalities to create a directed acyclic graph of GO terms, filter terms, and filter annotations.
        """
        if file_resources is None:
            file_resources = {
                "go-basic.obo": "http://purl.obolibrary.org/obo/go/go-basic.obo",
                "goa_human.gaf": "goa_human.gaf.gz",
                "goa_human_rna.gaf": "goa_human_rna.gaf.gz",
                "goa_human_isoform.gaf": "goa_human_isoform.gaf.gz",
            }
        super().__init__(
            path,
            file_resources,
            col_rename=col_rename,
            npartitions=npartitions,
            verbose=verbose,
        )

    def info(self):
        print("network {}".format(nx.info(self.network)))

    def load_dataframe(self, file_resources: Dict[str, TextIOWrapper], npartitions=None):
        # Annotations for each GO term
        go_annotations = pd.DataFrame.from_dict(dict(self.network.nodes(data=True)), orient='index')
        go_annotations["def"] = go_annotations["def"].apply(lambda x: x.split('"')[1] if isinstance(x, str) else None)
        go_annotations.index.name = "go_id"

        # Handle .gaf annotation files
        gaf_annotation_dfs = []
        for file in file_resources:
            if ".gaf" in file:
                go_lines = []
                for line in tqdm.tqdm(gafiterator(file_resources[file]), desc=file):
                    go_lines.append(line)
                gaf_annotation_dfs.append(pd.DataFrame(go_lines))

        if len(gaf_annotation_dfs):
            self.annotations: pd.DataFrame = pd.concat(gaf_annotation_dfs).reset_index(drop=True)
            self.annotations = self.annotations.rename(columns=self.COLUMNS_RENAME_DICT)

            self.annotations["Date"] = pd.to_datetime(self.annotations["Date"], )

            print("gaf_annotations:", self.annotations.columns.tolist())

        return go_annotations

    def load_network(self, file_resources):
        for file in file_resources:
            if ".obo" in file:
                network: nx.MultiDiGraph = obonet.read_obo(file_resources[file])
                network = network.reverse(copy=True)
                node_list = np.array(network.nodes)

        return network, node_list

    def annotation_train_val_test_split(self, train_date: str = "2017-06-15",
                                        valid_date: str = "2017-11-15",
                                        test_date: str = "2021-12-31",
                                        filter_evidence: List = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC'],
                                        groupby: List[str] = ["gene_name", "Qualifier"],
                                        filter_go_id: List[str] = None,
                                        ) -> Tuple[DataFrame, DataFrame, DataFrame]:
        annotations = self.annotations[self.annotations["Evidence"].isin(filter_evidence)]
        if filter_go_id is not None:
            annotations = annotations[annotations["go_id"].isin(filter_go_id)]

        # Split train/valid/test annotations
        train_anns = annotations[annotations["Date"] <= pd.to_datetime(train_date)]
        valid_anns = annotations[annotations["Date"] <= pd.to_datetime(valid_date)]

        test_anns = annotations.drop(index=valid_anns.index)
        valid_anns = valid_anns.drop(index=train_anns.index)
        if test_date:
            test_anns = test_anns[test_anns["Date"] <= pd.to_datetime(test_date)]

        # Set the source column (i.e. protein_id or gene_name), to be the first in groupby
        src_node_col = groupby[0]
        _exclude_single = lambda li: len(li) == 1 and "GO:0005515" in li

        outputs = []
        for anns in [train_anns, valid_anns, test_anns]:
            is_neg_ann = anns["Qualifier"].map(lambda li: "NOT" in li)
            # Convert `Qualifiers` entries of lists to strings
            anns = anns.apply({"Qualifier": lambda li: "".join([i for i in li if i != "NOT"]),
                               src_node_col: lambda x: x,
                               "go_id": lambda x: x,
                               "Evidence": lambda x: x})

            # Positive gene-GO annotations
            gene_anns: DataFrame = anns[~is_neg_ann].groupby(groupby).agg(go_id=("go_id", "unique"))

            # Negative gene-GO annotations
            neg_anns = anns[is_neg_ann].groupby(groupby).agg(neg_go_id=("go_id", "unique"))
            gene_anns["neg_go_id"] = neg_anns["neg_go_id"]
            gene_anns.drop(index=[""], inplace=True, errors="ignore")

            # Remove "GO:0005515" (protein binding) annotations for a gene if it's the gene's only annotation
            gene_anns.loc[gene_anns["go_id"].map(_exclude_single), "go_id"] = None
            gene_anns.drop(index=gene_anns.index[gene_anns.isna().all(1)], inplace=True)

            outputs.append(gene_anns)

        return tuple(outputs)


class InterPro(GeneOntology):
    def __init__(self, path="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/",
                 file_resources=None, col_rename=None, verbose=False, npartitions=None):
        """
        Args:
            path:
            file_resources:
            col_rename:
            verbose:
            npartitions:
        """

        if file_resources is None:
            file_resources = {}
            file_resources["entry.list"] = os.path.join(path, "entry.list")
            file_resources["protein2ipr.dat"] = os.path.join(path, "protein2ipr.dat.gz")
            file_resources["interpro2go"] = os.path.join(path, "interpro2go")
            file_resources["ParentChildTreeFile.txt"] = os.path.join(path, "ParentChildTreeFile.txt")

        super().__init__(path, file_resources, col_rename, verbose=verbose, npartitions=npartitions)

    def load_dataframe(self, file_resources: Dict[str, TextIOWrapper], npartitions=None):
        ipr_entries = pd.read_table(file_resources["entry.list"])
        ipr2go = self.parse_interpro2go(file_resources["interpro2go"])

        ipr_entries = ipr_entries.join(ipr2go.groupby('ENTRY_AC')["go_id"].unique(), on="ENTRY_AC")

        if "ppi_mat.npz" in file_resources:
            def get_pid_list(pid_list_file):
                with open(pid_list_file) as fp:
                    return [line.split()[0] for line in fp]

            net_pid_list = get_pid_list(file_resources["ppi_pid_list.txt"])
            node_feats = ssp.load_npz(file_resources["ppi_mat.npz"])
            df_feats = pd.DataFrame.sparse.from_spmatrix(node_feats, index=net_pid_list)
            self.annotations = df_feats
        else:
            self.annotations: dd.DataFrame = dd.read_table(
                file_resources["protein2ipr.dat"],
                names=['UniProtKB-AC', 'ENTRY_AC', 'ENTRY_NAME', 'accession', 'start', 'stop'],
                usecols=['UniProtKB-AC', 'ENTRY_AC'],
                dtype={'UniProtKB-AC': 'category', 'ENTRY_AC': 'category'})

        return ipr_entries

    def load_network(self, file_resources):
        for file in file_resources:
            if 'ParentChildTreeFile' in file and isinstance(file_resources[file], str) \
                and os.path.exists(file_resources[file]):
                network: nx.MultiDiGraph = self.parse_ipr_treefile(file_resources[file])
                node_list = np.array(network.nodes)

                return network, node_list
        return None, None

    def parse_ipr_treefile(self, lines: Union[Iterable[str], StringIO]) -> nx.MultiDiGraph:
        """Parse the InterPro Tree from the given file.
        Args:
            lines: A readable file or file-like
        """
        if isinstance(lines, str):
            lines = open(lines, 'r')

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

        return graph

    def parse_interpro2go(self, file: StringIO) -> DataFrame:
        if isinstance(file, str):
            file = open(file, 'r')

        def _process_line(line: str) -> Tuple[str, str, str]:
            pos = line.find('> GO')
            interpro_terms, go_term = line[:pos], line[pos:]
            interpro_id, interpro_name = interpro_terms.strip().split(' ', 1)
            go_name, go_id = go_term.split(';')
            go_desc = go_name.strip('> GO:')

            return (interpro_id.strip().split(':')[1], go_id.strip(), go_desc)

        tuples = [_process_line(line.strip()) for line in file if line[0] != '!']
        return pd.DataFrame(tuples, columns=['ENTRY_AC', "go_id", "go_desc"])


class UniProtGOA(GeneOntology):
    """Loads the GeneOntology database from https://www.ebi.ac.uk/GOA/ .

    Default path: "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/" .
    Default file_resources: {
        "goa_uniprot_all.gpi": "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpi.gz",
        "goa_uniprot_all.gaf": "goa_uniprot_all.gaf.gz",
    }
    """
    COLUMNS_RENAME_DICT = {
        "DB_Object_ID": "UniProtKB-AC",
        "DB_Object_Symbol": "gene_name",
        "GO_ID": "go_id",
    }

    def __init__(
        self,
        path="ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/",
        species="HUMAN",
        file_resources=None,
        col_rename=COLUMNS_RENAME_DICT,
        npartitions=0,
        verbose=False,
    ):
        """
        Loads the GeneOntology database from http://geneontology.org .

            Default path: "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/" .
            Default file_resources: {
                "goa_uniprot_all.gpi": "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpi.gz",
                "goa_uniprot_all.gaf": "goa_uniprot_all.gaf.gz",
            }

        Handles downloading the latest Gene Ontology obo and annotation data, preprocesses them. It provide
        functionalities to create a directed acyclic graph of GO terms, filter terms, and filter annotations.
        """
        self.species = species

        if file_resources is None:
            file_resources = {
                "go.obo": "http://current.geneontology.org/ontology/go.obo",
                "goa_human.gpi": os.path.join(species, "goa_human.gpi.gz"),
                "goa_human.gaf": os.path.join(species, "goa_human.gaf.gz"),
                "goa_human_isoform.gaf": os.path.join(species, "goa_human_isoform.gaf.gz"),
                # "goa_human_complex.gaf": os.path.join(species, "goa_human_complex.gaf.gz"),
            }
        super().__init__(
            path,
            file_resources,
            col_rename=col_rename,
            npartitions=npartitions,
            verbose=verbose,
        )


def gafiterator(handle):
    inline = handle.readline()
    if inline.strip().startswith("!gaf-version: 2"):
        # sys.stderr.write("gaf 2.0\n")
        return _gaf20iterator(handle)
    elif inline.strip() == "!gaf-version: 1.0":
        # sys.stderr.write("gaf 1.0\n")
        return _gaf10iterator(handle)
    else:
        return _gaf20iterator(handle)


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

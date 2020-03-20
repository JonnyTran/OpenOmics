import networkx as nx
import numpy as np
import obonet
import pandas as pd
from Bio.UniProt import GOA

from .base import Dataset


class GeneOntology(Dataset):
    COLUMNS_RENAME_DICT = {
        "DB_Object_Symbol": "gene_name",
        "DB_Object_ID": "gene_id",
        "GO_ID": "go_id"
    }

    def __init__(self, path="http://geneontology.org/gene-associations/",
                 file_resources=None, col_rename=COLUMNS_RENAME_DICT, npartitions=0):
        if file_resources is None:
            file_resources = {
                "go-basic.obo": "http://purl.obolibrary.org/obo/go/go-basic.obo",
                "goa_human.gaf": "goa_human.gaf.gz",
                "goa_human_rna.gaf": "goa_human_rna.gaf.gz",
                "goa_human_isoform.gaf": "goa_human_isoform.gaf.gz"
            }
        super(GeneOntology, self).__init__(path, file_resources, col_rename=col_rename, npartitions=npartitions)

        print("network {}".format(nx.info(self.network)))

    def load_dataframe(self, file_resources):
        go_annotation_dfs = []
        for file in file_resources:
            if ".gaf" in file:
                go_lines = []
                for line in GOA.gafiterator(file_resources[file]):
                    go_lines.append(line)
                go_annotation_dfs.append(pd.DataFrame(go_lines))

        go_annotations = pd.concat(go_annotation_dfs)

        for file in file_resources:
            if ".obo" in file:
                self.network = obonet.read_obo(file_resources[file])
                go_terms = pd.DataFrame.from_dict(self.network.nodes, orient='index', dtype="object")

                go_annotations["go_name"] = go_annotations["GO_ID"].map(go_terms["name"])
                go_annotations["namespace"] = go_annotations["GO_ID"].map(go_terms["namespace"])
                go_annotations["is_a"] = go_annotations["GO_ID"].map(go_terms["is_a"])

        return go_annotations

    def filter_terms(self, annotation, namespace=None):
        if namespace:
            biological_process_terms = set(
                self.df[self.df["namespace"] == "biological_process"]["go_id"].unique())
            filtered_annotation = annotation.str.split("|").map(
                lambda x: list(set(term for term in x if x in self.network) & biological_process_terms) if x else None)
        else:
            filtered_annotation = annotation.str.split("|").map(
                lambda x: [term for term in x if x in self.network] if x else None)

        return filtered_annotation

    def get_predecessor_terms(self, annotation: pd.Series):
        go_terms_parents = annotation.map(
            lambda terms: list({parent for term in terms for parent in self.network.predecessors(term)}) \
                if isinstance(terms, list) else None)
        return go_terms_parents

    def add_predecessor_terms(self, annotation: pd.Series, return_str=True):
        if annotation.dtypes == np.object and annotation.str.contains("\||;", regex=True).any():
            go_terms_annotations = annotation.str.split("|")
        else:
            go_terms_annotations = annotation

        go_terms_parents = go_terms_annotations + self.get_predecessor_terms(annotation)

        if return_str:
            go_terms_parents = go_terms_parents.map(
                lambda x: "|".join(x) if isinstance(x, list) else None)

        return go_terms_parents

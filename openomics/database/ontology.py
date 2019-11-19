
import os

import pandas as pd
from Bio.UniProt import GOA

from .annotation import Dataset


class GeneOntology(Dataset):
    def __init__(self, path, file_resources=None, col_rename=None, npartitions=0):

        if file_resources is None:
            file_resources = {}
            file_resources["gene_with_protein_product.txt"] = os.path.join(path, "gene_with_protein_product.txt")

        super(Dataset, self).__init__(path=path, file_resources=file_resources, col_rename=col_rename,
                                      npartitions=npartitions)

    def load_dataframe(self, file_resources):
        pass

    def get_GO_genes_info(self):
        lines = []
        with open(self.gene_ontology_file_path) as file:
            l = GOA.gafiterator(file)
            for line in l:
                lines.append(line)
        go_df = pd.DataFrame(lines)
        return go_df

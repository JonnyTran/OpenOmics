
import os

from openomics.database.base import Dataset


class GeneOntology(Dataset):
    def __init__(self, path, file_resources=None, col_rename=None, npartitions=0):

        if file_resources is None:
            file_resources = {}
            file_resources["gene_with_protein_product.txt"] = os.path.join(path, "gene_with_protein_product.txt")

        super().__init__(path, file_resources, col_rename=col_rename, npartitions=npartitions)

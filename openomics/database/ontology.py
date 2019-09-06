
import os

from openomics.database.annotation import Dataset


class GeneOntology(Dataset):
    def __init__(self, import_folder, file_resources=None, col_rename=None, npartitions=0):

        if file_resources is None:
            file_resources = {}
            file_resources["gene_with_protein_product.txt"] = os.path.join(import_folder, "gene_with_protein_product.txt")

        super().__init__(import_folder, file_resources, col_rename=col_rename, npartitions=npartitions)

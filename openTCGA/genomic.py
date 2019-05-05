import os

from openTCGA.expression import ExpressionData


class SomaticMutation(ExpressionData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "somaticMutation_geneLevel.txt")
        super().__init__(cancer_type, file_path)


class DNAMethylation(ExpressionData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "methylation_450.txt")
        super().__init__(cancer_type, file_path)


class CopyNumberVariation(ExpressionData):
    def __init__(self, cancer_type, folder_path):
        file_path = os.path.join(folder_path, "copyNumber.txt")
        super().__init__(cancer_type, file_path)
from openomics.database.annotation import Annotatable
from openomics.transcriptomics import ExpressionData


class SomaticMutation(ExpressionData, Annotatable):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super().__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                         transposed=transposed, log2_transform=log2_transform)

    def name(self):
        return __class__.__name__


class DNAMethylation(ExpressionData):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super().__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                         transposed=transposed, log2_transform=log2_transform)

    def name(self):
        return __class__.__name__


class CopyNumberVariation(ExpressionData, Annotatable):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super().__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                         transposed=transposed, log2_transform=log2_transform)

    def name(self):
        return __class__.__name__

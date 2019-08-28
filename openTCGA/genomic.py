from openTCGA.expression import ExpressionData

class SomaticMutation(ExpressionData):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super().__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                         transposed=transposed, log2_transform=log2_transform)


class DNAMethylation(ExpressionData):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super().__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                         transposed=transposed, log2_transform=log2_transform)


class CopyNumberVariation(ExpressionData):
    def __init__(self, cohort_name, index, file_path, columns, genes_col_name, transposed=True, log2_transform=False):
        super().__init__(cohort_name, index, file_path, columns=columns, genes_col_name=genes_col_name,
                         transposed=transposed, log2_transform=log2_transform)

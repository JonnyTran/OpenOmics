from openTCGA.expression import ExpressionData

class SomaticMutation(ExpressionData):
    def __init__(self, cohort_name, level, file_path, columns, index, transposed, log2_transform):
        super().__init__(cohort_name, level, file_path, columns=columns, index=index, transposed=transposed,log2_transform=log2_transform)


class DNAMethylation(ExpressionData):
    def __init__(self, cohort_name, level, file_path, columns, index, transposed, log2_transform):
        super().__init__(cohort_name, level, file_path, columns=columns, index=index, transposed=transposed, log2_transform=log2_transform)


class CopyNumberVariation(ExpressionData):
    def __init__(self, cohort_name, level, file_path, columns, index, transposed, log2_transform):
        super().__init__(cohort_name, level, file_path, columns=columns, index=index, transposed=transposed, log2_transform=log2_transform)

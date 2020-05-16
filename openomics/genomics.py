from openomics import ExpressionData
from openomics.database import Annotatable


class SingleNucleotideVariants(ExpressionData, Annotatable):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None,
                 sample_index_by="sample_index", log2_transform=False, npartitions=None):
        super(SingleNucleotideVariants, self).__init__(cohort_name, data, transposed=transposed, columns=columns,
                                                       gene_index_by=gene_index_by,
                                                       sample_index_by=sample_index_by, log2_transform=log2_transform,
                                                       npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__


class SomaticMutation(ExpressionData, Annotatable):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None,
                 sample_index_by="sample_index", log2_transform=False, npartitions=None):
        super(SomaticMutation, self).__init__(cohort_name, data, transposed=transposed, columns=columns,
                                              gene_index_by=gene_index_by,
                                              sample_index_by=sample_index_by, log2_transform=log2_transform,
                                              npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__


class DNAMethylation(ExpressionData, Annotatable):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None,
                 sample_index_by="sample_index", log2_transform=False, npartitions=None):
        super(DNAMethylation, self).__init__(cohort_name, data, transposed=transposed, columns=columns,
                                             gene_index_by=gene_index_by,
                                             sample_index_by=sample_index_by, log2_transform=log2_transform,
                                             npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__


class CopyNumberVariation(ExpressionData, Annotatable):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None,
                 sample_index_by="sample_index", log2_transform=False, npartitions=None):
        super(CopyNumberVariation, self).__init__(cohort_name, data, transposed=transposed, columns=columns,
                                                  gene_index_by=gene_index_by,
                                                  sample_index_by=sample_index_by, log2_transform=log2_transform,
                                                  npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__

from openomics import ExpressionData
from openomics.database import Annotatable


class SingleNucleotideVariants(ExpressionData, Annotatable):
    def __init__(
        self,
        data,
        transpose,
        gene_index,
        usecols=None,
        gene_level=None,
        sample_level="sample_index",
        transform_fn=None,
        dropna=False,
        npartitions=None,
        cohort_name=None,
    ):
        """
        Args:
            data:
            transpose:
            gene_index:
            usecols:
            gene_level:
            sample_level:
            transform_fn:
            dropna:
            npartitions:
            cohort_name:
        """
        super(SingleNucleotideVariants, self).__init__(
            data=data,
            transpose=transpose,
            gene_index=gene_index,
            usecols=usecols,
            gene_level=gene_level,
            sample_level=sample_level,
            transform_fn=transform_fn,
            dropna=dropna,
            npartitions=npartitions,
            cohort_name=cohort_name,
        )

    @classmethod
    def name(cls):
        return cls.__name__


class SomaticMutation(ExpressionData, Annotatable):
    def __init__(
        self,
        data,
        transpose,
        gene_index,
        usecols=None,
        gene_level=None,
        sample_level="sample_index",
        transform_fn=None,
        dropna=False,
        npartitions=None,
        cohort_name=None,
    ):
        """
        Args:
            data:
            transpose:
            gene_index:
            usecols:
            gene_level:
            sample_level:
            transform_fn:
            dropna:
            npartitions:
            cohort_name:
        """
        super(SomaticMutation, self).__init__(
            data=data,
            transpose=transpose,
            gene_index=gene_index,
            usecols=usecols,
            gene_level=gene_level,
            sample_level=sample_level,
            transform_fn=transform_fn,
            dropna=dropna,
            npartitions=npartitions,
            cohort_name=cohort_name,
        )

    @classmethod
    def name(cls):
        return cls.__name__


class DNAMethylation(ExpressionData, Annotatable):
    def __init__(
        self,
        data,
        transpose,
        gene_index,
        usecols=None,
        gene_level=None,
        sample_level="sample_index",
        transform_fn=None,
        dropna=False,
        npartitions=None,
        cohort_name=None,
    ):
        """
        Args:
            data:
            transpose:
            gene_index:
            usecols:
            gene_level:
            sample_level:
            transform_fn:
            dropna:
            npartitions:
            cohort_name:
        """
        super(DNAMethylation, self).__init__(
            data=data,
            transpose=transpose,
            gene_index=gene_index,
            usecols=usecols,
            gene_level=gene_level,
            sample_level=sample_level,
            transform_fn=transform_fn,
            dropna=dropna,
            npartitions=npartitions,
            cohort_name=cohort_name,
        )

    @classmethod
    def name(cls):
        return cls.__name__


class CopyNumberVariation(ExpressionData, Annotatable):
    def __init__(
        self,
        data,
        transpose,
        gene_index,
        usecols=None,
        gene_level=None,
        sample_level="sample_index",
        transform_fn=None,
        dropna=False,
        npartitions=None,
        cohort_name=None,
    ):
        """
        Args:
            data:
            transpose:
            gene_index:
            usecols:
            gene_level:
            sample_level:
            transform_fn:
            dropna:
            npartitions:
            cohort_name:
        """
        super(CopyNumberVariation, self).__init__(
            data=data,
            transpose=transpose,
            gene_index=gene_index,
            usecols=usecols,
            gene_level=gene_level,
            sample_level=sample_level,
            transform_fn=transform_fn,
            dropna=dropna,
            npartitions=npartitions,
            cohort_name=cohort_name,
        )

    @classmethod
    def name(cls):
        return cls.__name__

from .database.base import Annotatable
from .transcriptomics import Expression

__all__ = ['Protein']

class Protein(Expression, Annotatable):
    def __init__(self, data, transpose, gene_index=None, usecols=None, gene_level=None, sample_level="sample_index",
                 transform_fn=None, dropna=False, npartitions=None, cohort_name=None):
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
        super().__init__(data=data, transpose=transpose, gene_index=gene_index, usecols=usecols,
                         gene_level=gene_level, sample_level=sample_level, transform_fn=transform_fn,
                         dropna=dropna, npartitions=npartitions, cohort_name=cohort_name)

    @classmethod
    def name(cls):
        return cls.__name__


import io
import os
from glob import glob

import dask.dataframe as dd
import numpy as np
import pandas as pd
# from Bio.UniProt import GOA
from dask import delayed

from .database import Annotatable
from .utils.df import drop_duplicate_columns


class ExpressionData(object):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None, sample_index_by="sample_index",
                 log2_transform=False, dropna=False, npartitions=None):
        """
        .. class:: ExpressionData
        An abstract class that handles importing of any quantitative -omics data that is in a table format (e.g. csv, tsv, excel). Pandas will load the DataFrame from file with the user-specified columns and genes column name, then tranpose it such that the rows are samples and columns are gene/transcript/peptides.
        The user will also specify the index argument, which specifies if the genes are ensembl genes ID or gene name, or transcripts id/names. The user should be careful about choosing the right genes index which makes it easier to annotate functional, sequence, and interaction data to it.
        The dataframe should only contain numeric values besides the genes_col_name and the sample barcode id indices.
        Args:
            dropna:
            cohort_name (str): the unique cohort code name string
            data (str, byte-like, pandas.DataFrame):
                Path or file stream of the table file to import. If a pandas DataFrame is passed, then import this dataframe and skip preprocessing steps.
            transposed (bool): default True
                True if sample names are columns and rows are genes. False if the table has samples for row index, and gene names as columns.
            columns (str): a regex string
                A regex string to import column names from the table. Columns names imported are string match, separated by "|".
            gene_index_by (str): {"gene_id", "transcript_id", "peptide_id", "gene_name", "trascript_name", "peptide_name"}
                Chooses the level of the gene/transcript/peptide of the genes list in this expression data. The expression DataFrame's index will be renamed to this.
            sample_index_by (str): {"sample_index", "patient_index"}
                Chooses the level of the patient/sample/aliquot indexing.
            log2_transform (bool): default False
                Whether to log2 transform the expression values
            npartitions (int): [0-n], default 0
                If 0, then uses a Pandas DataFrame, if >1, then creates an off-memory Dask DataFrame with n partitions
        """
        self.cohort_name = cohort_name

        df = self.load_dataframe(data, transposed=transposed, columns=columns, genes_index=gene_index_by)
        self.expressions = self.preprocess_table(df, columns=columns, genes_index=gene_index_by, transposed=transposed,
                                                 dropna=dropna)
        if npartitions:
            self.expressions = dd.from_pandas(self.expressions, npartitions=npartitions)

        self.gene_index = gene_index_by
        self.sample_index = sample_index_by
        self.expressions.index.name = self.sample_index

        if log2_transform:
            self.expressions = self.expressions.applymap(self.log2_transform)

    def load_dataframe(self, data, transposed, columns, genes_index):
        """
        Reading table data inputs to create a DataFrame.

        Args:
            data: either a file path, a glob file path (e.g. "table-*.tsv"), a pandas.DataFrame, or a dask DataFrame.
            transposed: True if table oriented with samples
            columns:
            genes_index:

        Returns:

        """
        if isinstance(data, pd.DataFrame):
            df = data
            df = df.reset_index()
        elif isinstance(data, dd.DataFrame):
            df = data
        elif "*" in data:
            df = self.load_dataframe_glob(data, columns, genes_index, transposed)
        elif isinstance(data, io.StringIO):
            # TODO implement handling for multiple file ByteIO
            data.seek(0)  # Needed since the file was previous read to extract columns information
            df = pd.read_table(data)
        elif type(data) == str and os.path.isfile(data):
            df = pd.read_table(data, sep=None)
        else:
            raise IOError(data)
        return df

    def load_dataframe_glob(self, glob_path, columns, genes_index, transposed):
        # type: (str, str, str, bool) -> dd.DataFrame
        """

        :param glob_path:
        :param columns:
        :param genes_index:
        :param transposed:
        :return:
        """
        lazy_dataframes = []
        for file_path in glob(glob_path):
            df = delayed(pd.read_table)(file_path, )
            df = delayed(self.preprocess_table)(df, columns, genes_index, transposed, True)
            lazy_dataframes.append(df)

        return dd.from_delayed(lazy_dataframes, divisions=None)

    def preprocess_table(self, df, columns=None, genes_index=None, transposed=True, sort_index=False, dropna=True):
        # type: (pd.DataFrame, str, str, bool) -> pd.DataFrame
        """
        This function preprocesses the expression table files where columns are samples and rows are gene/transcripts
        Args:
            df (DataFrame): A Dask or Pandas DataFrame
            columns (str): A regular expression string for the column names to fetch.
            genes_index (str): The column name containing the gene/transcript names or id's.
            transposed: Default True. Whether to transpose the dataframe so columns are genes (features) and rows are samples.
        Returns:
            dataframe: a processed Dask DataFrame
        """

        # Filter columns
        if columns is not None:
            if genes_index not in columns:
                columns = columns + "|" + genes_index
            df = df.filter(regex=columns)

        # Cut TCGA column names to sample barcode, discarding aliquot info
        # df = df.rename(columns=lambda x: x[:16] if ("TCGA" in x) else x)

        # Drop duplicate sample names
        df = drop_duplicate_columns(df)

        # Drop NA geneID rows
        if dropna:
            df.dropna(axis=0, inplace=True)

        # Remove entries with unknown geneID
        if genes_index is not None:
            df = df[df[genes_index] != '?']
            df.set_index(genes_index, inplace=True)

        # Needed for Dask Delayed
        if sort_index == True:
            df.sort_index(axis=0, ascending=True, inplace=True)

        # Select only numerical columns
        df = df.select_dtypes(include="number")

        # Transpose dataframe to sample rows and gene columns
        if transposed:
            df = df.T

        # Drop duplicate genes
        df = drop_duplicate_columns(df)

        return df

    def set_genes_index(self, index: str, old_index: str):
        assert isinstance(self, Annotatable) and isinstance(self, ExpressionData)
        # Change gene name columns in expressions
        rename_dict = self.get_rename_dict(from_index=old_index, to_index=index)
        self.expressions.rename(columns=rename_dict, inplace=True)
        self.gene_index = index

        # Change index name in annotation
        self.set_index(index)

    def log2_transform(self, x):
        return np.log2(x + 1)

    def drop_genes(self, genes_to_drop):
        self.expressions.drop(genes_to_drop, axis=1, inplace=True)
        for gene in genes_to_drop:
            self.features.remove(gene)

    @classmethod
    def name(cls):
        raise NotImplementedError

    def get_genes_list(self, level=None):
        index = self.expressions.columns

        if isinstance(index, pd.MultiIndex):
            return index.get_level_values(self.gene_index if level is None else level)
        else:
            return index

    def get_samples_list(self, level=None):
        index = self.expressions.index
        if isinstance(index, pd.MultiIndex):
            return index.get_level_values(self.gene_index if level is None else level)
        else:
            return index
        return index

    samples = property(get_samples_list)
    features = property(get_genes_list)


class LncRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None, sample_index_by="sample_index",
                 log2_transform=False, dropna=False, npartitions=None):
        super(LncRNA, self).__init__(cohort_name, data=data, transposed=transposed, columns=columns,
                                     gene_index_by=gene_index_by, sample_index_by=sample_index_by,
                                     log2_transform=log2_transform, dropna=dropna, npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__


class MessengerRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None, sample_index_by="sample_index",
                 log2_transform=False, dropna=False, npartitions=None):
        super(MessengerRNA, self).__init__(cohort_name, data=data, transposed=transposed, columns=columns,
                                           gene_index_by=gene_index_by, sample_index_by=sample_index_by,
                                           log2_transform=log2_transform, dropna=dropna, npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__


class MicroRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, data, transposed, columns=None, gene_index_by=None, sample_index_by="sample_index",
                 log2_transform=False, dropna=False, npartitions=None):
        super(MicroRNA, self).__init__(cohort_name, data=data, transposed=transposed, columns=columns,
                                       gene_index_by=gene_index_by, sample_index_by=sample_index_by,
                                       log2_transform=log2_transform, dropna=dropna, npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__

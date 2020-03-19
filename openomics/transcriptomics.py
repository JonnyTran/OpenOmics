import io
import os
from glob import glob

import dask.dataframe as dd
import numpy as np
import pandas as pd
# from Bio.UniProt import GOA
from dask import delayed

from .database import Annotatable


class ExpressionData(object):
    def __init__(self, cohort_name, file_path, columns=None, genes_col_name=None, gene_index_by=None,
                 sample_index_by="sample_index", transposed=True, log2_transform=False, npartitions=None):
        """
        .. class:: ExpressionData
        An abstract class that handles importing of any quantitative -omics data that is in a table format (e.g. csv, tsv, excel). Pandas will load the DataFrame from file with the user-specified columns and genes column name, then tranpose it such that the rows are samples and columns are gene/transcript/peptides.
        The user will also specify the index argument, which specifies if the genes are ensembl genes ID or gene name, or transcripts id/names. The user should be careful about choosing the right genes index which makes it easier to annotate functional, sequence, and interaction data to it.
        The dataframe should only contain numeric values besides the genes_col_name and the sample barcode id indices.
        Args:
            cohort_name (str): the unique cohort code name string
            file_path (str, byte-like, Pandas DataFrame):
                Path or file stream of the table file to import. If a pandas DataFrame is passed, then import this dataframe and skip preprocessing steps.
            columns (str): a regex string
                A regex string to import column names from the table. Columns names imported are string match, separated by "|".
            genes_col_name (str):
                Index column name of expression data which is used to index the genes list
            gene_index_by (str): {"gene_id", "transcript_id", "peptide_id", "gene_name", "trascript_name", "peptide_name"}
                Chooses the level of the gene/transcript/peptide of the genes list in this expression data. The expression DataFrame's index will be renamed to this.
            sample_index_by (str): {"sample_index", "patient_index"}
                Chooses the level of the patient/sample/aliquot indexing.
            transposed (bool): default True
                True if sample names are columns and rows are genes. False if the table has samples for row index, and gene names as columns.
            log2_transform (bool): default False
                Whether to log2 transform the expression values
            npartitions (int): [0-n], default 0
                If 0, then uses a Pandas DataFrame, if >1, then creates an off-memory Dask DataFrame with n partitions
        """
        self.cohort_name = cohort_name

        if isinstance(file_path, pd.DataFrame):
            df = file_path
        elif "*" in file_path:
            df = self.preprocess_table_glob(file_path, columns, genes_col_name, transposed)
        elif isinstance(file_path, io.StringIO):
            # TODO implement handling for multiple file ByteIO
            file_path.seek(0)  # Needed since the file was previous read to extract columns information
            df = pd.read_table(file_path)
        elif type(file_path) == str and os.path.isfile(file_path):
            df = pd.read_table(file_path, sep=None)
        else:
            raise IOError(file_path)

        self.expressions = self.preprocess_table(df, columns, genes_col_name, transposed)
        if npartitions:
            self.expressions = dd.from_pandas(self.expressions, npartitions=npartitions)

        self.gene_index = gene_index_by
        self.sample_index = sample_index_by
        self.expressions.index.name = self.sample_index

        if log2_transform:
            self.expressions = self.expressions.applymap(self.log2_transform)

    def set_genes_index(self, index, old_index):
        assert isinstance(self, Annotatable) and isinstance(self, ExpressionData)
        # Change gene name columns in expressions
        rename_dict = self.get_rename_dict(from_index=old_index, to_index=index)
        self.expressions.rename(columns=rename_dict, inplace=True)
        self.gene_index = index

        # Change index name in annotation
        self.set_index(index)

    def preprocess_table(self, df, columns=None, genes_index=None, transposed=True, sort_index=False):
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
        df = df.rename(columns=lambda x: x[:16] if ("TCGA" in x) else x)

        # Drop duplicate sample names
        _, i = np.unique(df.columns, return_index=True)
        df = df.iloc[:, i]

        # Drop NA geneID rows
        df.dropna(axis=0, inplace=True)

        # Remove entries with unknown geneID
        if genes_index is not None:
            df = df[df[genes_index] != '?']
            df.set_index(genes_index, inplace=True)

        # Needed for Dask Delayed
        if sort_index == True:
            df.sort_index(axis=0, ascending=True, inplace=True)

        # Transpose dataframe to patient rows and geneID columns
        if transposed:
            df = df.T

        # Drop duplicate columns names (Gene symbols with same name)
        _, i = np.unique(df.columns, return_index=True)
        df = df.iloc[:, i]

        return df

    def preprocess_table_glob(self, glob_path, columns, genes_index, transposed):
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

    def log2_transform(self, x):
        return np.log2(x + 1)

    def drop_genes(self, genes_to_drop):
        self.expressions.drop(genes_to_drop, axis=1, inplace=True)
        for gene in genes_to_drop:
            self.features.remove(gene)

    @classmethod
    def name(cls):
        raise NotImplementedError

    def get_genes_list(self):
        return self.expressions.columns

    def get_samples_list(self):
        return self.expressions.index

    samples = property(get_samples_list)
    features = property(get_genes_list)


class LncRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, file_path, columns, genes_col_name, gene_index_by=None,
                 sample_index_by="sample_barcode",
                 transposed=True, log2_transform=False, npartitions=0):
        super(LncRNA, self).__init__(cohort_name, file_path=file_path, columns=columns, genes_col_name=genes_col_name,
                                     gene_index_by=gene_index_by, sample_index_by=sample_index_by,
                                     transposed=transposed,
                                     log2_transform=log2_transform, npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__

    def preprocess_table(self, df, columns, genes_index, transposed):
        # type: (pd.DataFrame, str, str, bool) -> pd.DataFrame
        """
        Preprocess LNCRNA expression file obtained from TANRIC MDAnderson, and replace ENSEMBL gene ID to HUGO gene names (HGNC). This function overwrites the GenomicData.process_expression_table() function which processes TCGA-Assembler data.

        TANRIC LNCRNA expression values are log2 transformed
        :param transposed:
        """

        # Replacing ENSG Gene ID to the lncRNA gene symbol name
        df[genes_index] = df[genes_index].str.replace("[.].*", "")  # Removing .# ENGS gene version number at the end
        df = df[~df[genes_index].duplicated(keep='first')]  # Remove duplicate genes

        # Drop NA gene rows
        df.dropna(axis=0, inplace=True)

        # Transpose matrix to patients rows and genes columns
        df.index = df[genes_index]
        df = df.T.iloc[1:, :]

        # Change index string to bcr_sample_barcode standard
        def change_patient_barcode(s):
            if "Normal" in s:
                return s[s.find('TCGA'):] + "-11A"
            elif "Tumor" in s:
                return s[s.find('TCGA'):] + "-01A"
            else:
                return s

        df.index = df.index.map(change_patient_barcode)
        df.index.name = "gene_id"

        return df

    # def get_lncBase_miRNA_lncRNA_predicted_interactions_edgelist(self, data=True, rename_dict=None):
    #     records = []
    #     for record in SeqIO.parse(
    #             self.lncBase_predicted_interactions_file_path,
    #             "fasta"):
    #         records.append(record.id.split(","))
    #     lncbase_df = pd.DataFrame(records, columns=["Transcript_ID", "geneId", "mirna", "miTG-score"])
    #
    #     lncbase_df["miTG-score"] = lncbase_df["miTG-score"].astype(float)
    #     lncbase_df = lncbase_df[lncbase_df["miTG-score"] >= 0.9]
    #
    #     lncbase_df["geneId"] = lncbase_df["geneId"].str.replace("(\(.+|\))", "")
    #     lncbase_df["mirna"] = lncbase_df["mirna"].str.replace("(\(.+|\))", "")
    #     lncbase_df["mirna"] = lncbase_df["mirna"].str.lower()
    #     lncbase_df["mirna"] = lncbase_df["mirna"].str.replace("-3p.*|-5p.*", "")
    #
    #     lncBase_lncRNA_miRNA_network = nx.from_pandas_edgelist(lncbase_df, source='mirna', target='geneId',
    #                                                            edge_attr=["miTG-score"],
    #                                                            create_using=nx.DiGraph())
    #     if rename_dict is not None:
    #         lncBase_lncRNA_miRNA_network = nx.relabel_nodes(lncBase_lncRNA_miRNA_network, rename_dict)
    #
    #     return lncBase_lncRNA_miRNA_network.edges(data=data)



class MessengerRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, file_path, columns, genes_col_name, gene_index_by=None,
                 sample_index_by="sample_barcode",
                 transposed=True,
                 log2_transform=False, npartitions=0):
        super(MessengerRNA, self).__init__(cohort_name, file_path=file_path, columns=columns,
                                           genes_col_name=genes_col_name, gene_index_by=gene_index_by,
                                           sample_index_by=sample_index_by,
                                           transposed=transposed, log2_transform=log2_transform,
                                           npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__



class MicroRNA(ExpressionData, Annotatable):
    def __init__(self, cohort_name, file_path, columns=None, genes_col_name=None, gene_index_by=None,
                 sample_index_by="sample_barcode",
                 transposed=True,
                 log2_transform=False, npartitions=0):
        super(MicroRNA, self).__init__(cohort_name, file_path=file_path, columns=columns, genes_col_name=genes_col_name,
                                       gene_index_by=gene_index_by, sample_index_by=sample_index_by,
                                       transposed=transposed,
                                       log2_transform=log2_transform, npartitions=npartitions)

    @classmethod
    def name(cls):
        return cls.__name__


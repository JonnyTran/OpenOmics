import copy
import difflib
import gzip
import os
import zipfile
from abc import abstractmethod

import filetype
import pandas as pd
import rarfile
import validators

from openomics.utils.df import concat_uniques
from openomics.utils.io import get_pkg_data_filename


class Dataset(object):
    COLUMNS_RENAME_DICT = None  # Needs initialization since subclasses may use this field

    def __init__(self, path, file_resources=None, col_rename=None, npartitions=0, verbose=False):
        """
        This is an abstract class used to instantiate a database given a folder containing various file resources. When creating a Database class, the load_data function is called where the file resources are load as a DataFrame and performs necessary processings. This class provides an interface for RNA classes to annotate various genomic annotation, functional annotation, sequences, and disease associations.
        Args:
            path (str):
                The folder or url path containing the data file resources. If url path, the files will be downloaded and cached to the user's home folder (at ~/.astropy/).
            file_resources (dict): default None,
                Used to list required files for preprocessing of the database. A dictionary where keys are required filenames and value are file paths. If None, then the class constructor should automatically build the required file resources dict.
            col_rename (dict): default None,
                A dictionary to rename columns in the data table. If None, then automatically load defaults.
            npartitions (int): [0-n], default 0
                If 0, then uses a Pandas DataFrame, if >1, then creates an off-memory Dask DataFrame with n partitions
        """
        self.verbose = verbose

        self.validate_file_resources(file_resources, path)

        self.df = self.load_dataframe(file_resources)
        self.df = self.df.reset_index()
        if col_rename is not None:
            self.df = self.df.rename(columns=col_rename)

        self.info() if verbose else None

    def info(self):
        print("{}: {}".format(self.name(), self.df.columns.tolist()))

    def validate_file_resources(self, file_resources, path, verbose=False):
        if validators.url(path):
            for filename, filepath in copy.copy(file_resources).items():
                data_file = get_pkg_data_filename(path, filepath,
                                                  verbose=verbose)  # Download file and replace the file_resource path
                filetype_ext = filetype.guess(data_file)

                if filetype_ext is None:  # This if-clause is needed incase when filetype_ext is None, causing the next clause to fail
                    file_resources[filename] = data_file  # Returns the

                elif filetype_ext.extension == 'gz':
                    file_resources[filename] = gzip.open(data_file, 'rt')

                elif filetype_ext.extension == 'zip':
                    zf = zipfile.ZipFile(data_file, 'r')
                    for f in zf.infolist():
                        file_resources[filename] = zf.open(f.filename, mode='r')

                elif filetype_ext.extension == 'rar':
                    rf = rarfile.RarFile(data_file, 'r')
                    for f in rf.infolist():
                        file_resources[filename] = rf.open(f.filename, mode='r')
                else:
                    file_resources[filename] = data_file

        elif os.path.isdir(path) and os.path.exists(path):
            for _, filepath in file_resources.items():
                if not os.path.exists(filepath):
                    raise IOError(filepath)
        else:
            raise IOError(path)

        self.data_path = path
        self.file_resources = file_resources

    def close(self):
        # Close opened file resources
        for filename, filepath in self.file_resources.items():
            if type(self.file_resources[filename]) != str:
                self.file_resources[filename].close()

    @abstractmethod
    def load_dataframe(self, file_resources):
        # type: (dict) -> pd.DataFrame
        """
        Handles data preprocessing given the file_resources input, and returns a DataFrame.

        Args:
            file_resources (dict): A dict with keys as filenames and values as full file path.
            **kwargs: Optional
        """
        raise NotImplementedError

    @classmethod
    def name(cls):
        return cls.__name__

    def list_databases(self):
        return DEFAULT_LIBRARIES

    def get_annotations(self, index, columns):
        # type: (str, List[str]) -> Union[pd.DataFrame, dd.DataFrame]
        """
        Returns the Database's DataFrame such that it's indexed by :param index:, which then applies a groupby operation
        and aggregates all other columns by concatenating all unique values.

        operation aggregates
        Args:
            index (str): The index column name of the Dataframe
            columns (list): a list of column names

        Returns:
            df (DataFrame): A dataframe to be used for annotation

        """
        if columns is not None:
            if index in columns:
                df = self.df[columns]
                columns.pop(columns.index(index))
            else:
                df = self.df[columns + [index]]
        else:
            raise Exception(
                "The columns argument must be a list such that it's subset of the following columns in the dataframe",
                self.df.columns.tolist())

        if index != self.df.index.name and index in self.df.columns:
            df = df.set_index(index)

        # Groupby index, and Aggregate by all columns by concatenating unique values
        df = df.groupby(index).agg({k: concat_uniques for k in columns})

        if df.index.duplicated().sum() > 0:
            raise ValueError("DataFrame must not have duplicates in index")
        return df

    def get_expressions(self, index):
        return self.df.groupby(
            index).median()  # TODO if index by gene, aggregate medians of transcript-level expressions

    @abstractmethod
    def get_rename_dict(self, from_index, to_index):
        """
        Used to retrieve a lookup dictionary to convert from one index to another, e.g., gene_id to gene_name

        Args:
            from_index: an index on the DataFrame for key
            to_index: an index on the DataFrame for value

        Returns
            rename_dict (dict): a rename dict
        """
        raise NotImplementedError


# from .sequence import SequenceDataset
class Annotatable(object):
    """
    This class provides an interface for the omics to annotate external data downloaded from various databases. These data will be imported as attribute information to the genes, or interactions between the genes.
    """

    def __init__(self):
        pass

    def get_annotations(self):
        if hasattr(self, "annotations"):
            return self.annotations
        else:
            raise Exception("{} must run initialize_annotations() first.".format(self.name()))

    def get_annotation_expressions(self):
        if hasattr(self, "annotation_expressions"):
            return self.annotation_expressions
        else:
            raise Exception("{} must run annotate_expressions() first.".format(self.name()))

    def initialize_annotations(self, gene_list, index):
        if gene_list is None:
            gene_list = self.get_genes_list()

        self.annotations = pd.DataFrame(index=gene_list)
        self.annotations.index.name = index

    def annotate_genomics(self, database, index, columns, fuzzy_match=False):
        # type: (Dataset, str, List[str], bool) -> None
        """
        Performs a left outer join between the annotation and Database's DataFrame, on the index key. The index argument must be column present in both DataFrames.
        If there exists overlapping column in the join, then the fillna() is used to fill NaN values in the old column with non-NaN values from the new column.
        Args:
            database (openomics.annotation.Database): Database which contains annotation
            index (str): The column name which exists in both the annotation and Database's DataFrame
            columns (list): a list of column name to join to the annotation
            fuzzy_match (bool): default False. Whether to join the annotation by applying a fuzzy match on the index with difflib.get_close_matches(). It is very computationally expensive and thus should only be used sparingly.
        """
        database_df = database.get_annotations(index, columns)
        if fuzzy_match:
            database_df.index = database_df.index.map(lambda x: difflib.get_close_matches(x, self.annotations.index)[0])

        if index == self.annotations.index.name:
            self.annotations = self.annotations.join(database_df, on=index, rsuffix="_")
        else:
            if type(self.annotations.index) == pd.MultiIndex:
                old_index = self.annotations.index.names
            else:
                old_index = self.annotations.index.name

            self.annotations = self.annotations.reset_index()
            self.annotations.set_index(index, inplace=True)
            self.annotations = self.annotations.join(database_df, on=index, rsuffix="_").reset_index()
            self.annotations.set_index(old_index, inplace=True)

        # Merge columns if the database DataFrame has overlapping columns with existing column
        duplicate_columns = [col for col in self.annotations.columns if col[-1] == "_"]
        for new_col in duplicate_columns:
            old_col = new_col.strip("_")
            self.annotations[old_col].fillna(self.annotations[new_col], inplace=True, axis=0)
            self.annotations.drop(columns=new_col, inplace=True)


    def annotate_sequences(self, database, index, agg_sequences="longest", omic=None, **kwargs):
        # type: (Dataset, str, str) -> None
        # assert isinstance(database, SequenceDataset)
        if omic is None:
            omic = self.name()

        sequences_entries = database.get_sequences(index=index, omic=omic, agg_sequences=agg_sequences, **kwargs)

        if type(self.annotations.index) == pd.MultiIndex:
            self.annotations['Transcript sequence'] = self.annotations.index.get_level_values(index).map(
                sequences_entries)
        else:
            self.annotations["Transcript sequence"] = self.annotations.index.map(sequences_entries)

    def annotate_expressions(self, database, index, fuzzy_match=False):
        self.annotation_expressions = pd.DataFrame(index=self.annotations.index)

        if self.annotations.index.name == index:
            self.annotation_expressions = self.annotation_expressions.join(
                database.get_expressions(index=index))
        else:
            raise Exception("index argument must be one of", database.df.index)

    def annotate_interactions(self, database, index):
        # type: (Interactions, str) -> None
        raise NotImplementedError

    def annotate_diseases(self, database, index):
        # type: (DiseaseAssociation, str) -> None
        self.annotations["disease_associations"] = self.annotations.index.map(
            database.get_disease_assocs(index=index, ))

    def set_index(self, new_index):
        self.annotations[new_index].fillna(self.annotations.index.to_series(), axis=0, inplace=True)
        self.annotations = self.annotations.reset_index().set_index(new_index)

    def get_rename_dict(self, from_index, to_index):
        dataframe = self.annotations.reset_index()
        dataframe = dataframe[dataframe[to_index].notnull()]
        return pd.Series(dataframe[to_index].values,
                         index=dataframe[from_index]).to_dict()


DEFAULT_LIBRARIES = ["10KImmunomes"
                     "BioGRID"
                     "CCLE"
                     "DisGeNET"
                     "ENSEMBL"
                     "GENCODE"
                     "GeneMania"
                     "GeneOntology"
                     "GlobalBiobankEngine"
                     "GTEx"
                     "HMDD_miRNAdisease"
                     "HPRD_PPI"
                     "HUGO_Gene_names"
                     "HumanBodyMapLincRNAs"
                     "IntAct"
                     "lncBase"
                     "LNCipedia"
                     "LncReg"
                     "lncRInter"
                     "lncrna2target"
                     "lncRNA_data_repository"
                     "lncrnadisease"
                     "lncRNome"
                     "mirbase"
                     "miRTarBase"
                     "NHLBI_Exome_Sequencing_Project"
                     "NONCODE"
                     "NPInter"
                     "PIRD"
                     "RegNetwork"
                     "RISE_RNA_Interactions"
                     "RNAcentral"
                     "StarBase_v2.0"
                     "STRING_PPI"
                     "TargetScan"]

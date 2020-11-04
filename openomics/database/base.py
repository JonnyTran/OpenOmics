import copy
import difflib
import gzip
import os
from abc import abstractmethod

import filetype
import pandas as pd
import validators

from openomics.utils.df import concat_uniques
from openomics.utils.io import get_pkg_data_filename


class Dataset(object):
    COLUMNS_RENAME_DICT = None  # Needs initialization since subclasses may use this field

    def __init__(self, path, file_resources=None, col_rename=None, npartitions=0):
        """
        This is an abstract class used to instantiate a database given a folder containing various file resources. When creating a Database class, the load_data function is called where the file resources are load as a DataFrame and performs necessary processings. This class provides an interface for RNA classes to annotate various genomic annotations, functional annotations, sequences, and disease associations.
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
        # Download data from ftp/html path
        if validators.url(path) or any([validators.url(file) for file in file_resources]):
            for filename, filepath in copy.copy(file_resources).items():
                # Download the files and replace the file_resource paths
                data_file = get_pkg_data_filename(path, filepath)
                extension = filetype.guess(data_file).extension

                if extension == 'gz':
                    file_resources[filename] = gzip.open(data_file, 'r')
                else:
                    file_resources[filename] = data_file
            print(file_resources)

        # If using local data
        elif os.path.isdir(path) and os.path.exists(path):
            for _, filepath in file_resources.items():
                if not os.path.exists(filepath):
                    raise IOError(filepath)
        else:
            raise IOError(path)

        self.data_path = path
        self.file_resources = file_resources
        self.df = self.load_dataframe(file_resources)
        self.df = self.df.reset_index()
        if col_rename is not None:
            self.df = self.df.rename(columns=col_rename)
        print("{}: {}".format(self.name(), self.df.columns.tolist()))

        # Close opened file resources
        for filename, filepath in file_resources.items():
            if type(file_resources[filename]) != str:
                file_resources[filename].close()

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

    def get_sequences(self, index, omic=None):
        """
        Returns a dictionary where keys are
        Args:
            omic (str): {"lncRNA", "microRNA", "messengerRNA"}
            index (str): {"gene_id", "gene_name", "transcript_id", "transcript_name"}
                The index
        """
        raise NotImplementedError


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
            raise Exception("Must run initialize_annotations() first.")

    def initialize_annotations(self, gene_list, index):
        if gene_list is None:
            gene_list = self.get_genes_list()

        self.annotations = pd.DataFrame(index=gene_list)
        self.annotations.index.name = index

    def annotate_genomics(self, database, index, columns, fuzzy_match=False):
        # type: (Dataset, str, List[str], bool) -> None
        """
        Performs a left outer join between the annotations and Database's DataFrame, on the index key. The index argument must be column present in both DataFrames.
        If there exists overlapping column in the join, then the fillna() is used to fill NaN values in the old column with non-NaN values from the new column.
        Args:
            database (openomics.annotation.Database): Database which contains annotations
            index (str): The column name which exists in both the annotations and Database's DataFrame
            columns (list): a list of column name to join to the annotations
            fuzzy_match (bool): default False. Whether to join the annotations by applying a fuzzy match on the index with difflib.get_close_matches(). It is very computationally expensive and thus should only be used sparingly.
        """
        database_annotations = database.get_annotations(index, columns)
        if fuzzy_match:
            database_annotations.index = database_annotations.index.map(
                lambda x: difflib.get_close_matches(x, self.annotations.index)[0])

        if index == self.annotations.index.name:
            self.annotations = self.annotations.join(database_annotations, on=index, rsuffix="_")
        else:
            old_index = self.annotations.index.name
            self.annotations = self.annotations.reset_index()
            self.annotations.set_index(index, inplace=True)
            self.annotations = self.annotations.join(database_annotations, on=index, rsuffix="_")
            self.annotations = self.annotations.reset_index()
            self.annotations.set_index(old_index, inplace=True)

        # Merge columns if the database DataFrame has overlapping columns with existing column
        duplicate_columns = [col for col in self.annotations.columns if col[-1] == "_"]
        for col in duplicate_columns:
            self.annotations[col.strip("_")].fillna(self.annotations[col], inplace=True, axis=0)
            self.annotations.drop(columns=col, inplace=True)

    def annotate_sequences(self, database, index, omic):
        # type: (Dataset, str, str) -> None
        self.annotations["Transcript sequence"] = self.annotations.index.map(
            database.get_sequences(index=index, omic=omic))

    def annotate_interactions(self, database, index):
        # type: (Dataset, str) -> None
        raise NotImplementedError

    def annotate_diseases(self, database, index):
        # type: (Dataset, str) -> None
        raise NotImplementedError

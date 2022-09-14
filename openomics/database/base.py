import copy
import difflib
import logging
import os
import warnings
from abc import ABC, abstractmethod
from os.path import exists, join
from typing import Dict, Union, Any
from typing import List

import dask.dataframe as dd
import filetype
import pandas as pd
import validators

from openomics.utils.df import concat_uniques, concat
from openomics.utils.io import get_pkg_data_filename, decompress_file


class Database(object):
    """This is a base class used to instantiate an external Database given a a set
    of files from either local files or URLs. When creating a Database class, the
    `load_dataframe()` function is called where the file_resources are used to
    load (Pandas or Dask) DataFrames, then performs data wrangling to yield a
    dataframe at `self.data` . This class also provides an interface for -omics
    tables, e.g. `ExpressionData` , to annotate various annotations,
    expressions, sequences, and disease associations.
    """
    data: pd.DataFrame
    COLUMNS_RENAME_DICT = None  # Needs initialization since subclasses may use this field to rename columns in dataframes.

    def __init__(
        self,
        path,
        file_resources=None,
        col_rename=None,
        blocksize=None,
        verbose=False,
    ):
        """
        Args:
            path: The folder or url path containing the data file resources. If
                url path, the files will be downloaded and cached to the user's
                home folder (at ~/.astropy/).
            file_resources: Used to list required files for preprocessing of the
                database. A dictionary where keys are required filenames and
                value are file paths. If None, then the class constructor should
                automatically build the required file resources dict.
            col_rename (dict): default None, A dictionary to rename columns in
                the data table. If None, then automatically load defaults.
            blocksize (int): [0-n], default 0 If 0, then uses a Pandas
                DataFrame, if >1, then creates an off-memory Dask DataFrame with
                partitions where each partion contains `blocksize` rows.
            verbose (bool): Default False.
        """
        if blocksize:
            assert isinstance(blocksize, bool) or (isinstance(blocksize, int) and blocksize > 10), \
                f"blocksize ({blocksize}) is too small and will cause a huge overhead in Dask dataframes."

        self.blocksize = blocksize
        self.verbose = verbose
        self.data_path = path

        self.file_resources = self.load_file_resources(path, file_resources=file_resources, verbose=verbose)
        self.data = self.load_dataframe(self.file_resources, blocksize=blocksize)
        if col_rename is not None:
            self.data = self.data.rename(columns=col_rename)
            if self.data.index.name in col_rename:
                self.data.index.name = col_rename[self.data.index.name]

    def __repr__(self):
        return "{}: {}".format(self.name(), self.data.columns.tolist())

    def load_file_resources(self, base_path, file_resources, verbose=False) -> Dict[str, Any]:
        """For each file in file_resources, download the file if path+file is a
        URL or load from disk if a local path. Additionally unzip or unrar if
        the file is compressed.

        Args:
            base_path (str): The folder or url path containing the data file
                resources. If a url path, the files will be downloaded and
                cached to the user's home folder (at ~/.astropy/).
            file_resources (dict): default None, Used to list required files for
                preprocessing of the database. A dictionary where keys are
                required filenames and value are file paths. If None, then the
                class constructor should automatically build the required file
                resources dict.
            verbose:
        """
        file_resources_new = copy.copy(file_resources)

        for filename, filepath in file_resources.items():
            # Remote database file URL
            if validators.url(filepath) or validators.url(join(base_path, filepath)):
                filepath = get_pkg_data_filename(base_path, filepath)
                filepath_ext = filetype.guess(filepath)

            # Local database path
            elif exists(filepath) or exists(join(base_path, filepath)):
                if isinstance(filepath, str) and not exists(filepath):
                    if exists(os.path.join(base_path, filepath)):
                        filepath = os.path.join(base_path, filepath)
                    else:
                        warnings.warn(f"`base_path` is a local file directory, so all file_resources must be local. "
                                      f"Cannot use `filepath` = {filepath} with `base_path` = {base_path}")
                        continue

                filepath_ext = filetype.guess(filepath)

            else:
                # file_path is an external file outside of `base_path`
                filepath_ext = None

            # Update filepath on uncompressed file
            file_resources_new[filename] = filepath

            if filepath_ext:
                file, new_filename = decompress_file(filepath, filename, file_ext=filepath_ext)
                file_resources_new[new_filename] = file

        # if validators.url(base_path) or any(validators.url(filepath) for filepath in file_resources.values()):
        #     for filename, filepath in file_resources.items():
        #         # Download file (if not already cached) and replace the file_resource path
        #         if isinstance(filepath, str):
        #             filepath = get_pkg_data_filename(base_path, filepath)
        #             filepath_ext = filetype.guess(filepath)
        #         else:
        #             filepath_ext = None
        #
        #         uncomp_filename, file = decompress_file(filepath, filename, filepath_ext)
        #         file_resources_new[uncomp_filename] = file
        #
        # # Local database path
        # elif os.path.isdir(base_path) and exists(base_path):
        #     for filename, filepath in file_resources.items():
        #         if not exists(filepath):
        #             warnings.warn(f"`base_path` is a local file directory, so all file_resources must be local. "
        #                   f"`filepath` = {filepath}")
        #             continue
        #
        #         if isinstance(filepath, str):
        #             filepath_ext = filetype.guess(filepath)
        #         else:
        #             filepath_ext = None
        #
        #         uncomp_filename, file = decompress_file(filepath, filename, filepath_ext)
        #         file_resources_new[uncomp_filename] = file
        # else:
        #     raise IOError(
        #         f"`base_path` {base_path} not supported. Must be either a remote URL point or a local directory path.")

        logging.info(f"{self.name()} file_resources: {file_resources_new}")
        return file_resources_new


    def close(self):
        # Close opened file resources
        for filename, filepath in self.file_resources.items():
            if hasattr(self.file_resources[filename], 'close'):
                self.file_resources[filename].close()

    @abstractmethod
    def load_dataframe(self, file_resources: Dict[str, str], blocksize: int = None) -> pd.DataFrame:
        """Handles data preprocessing given the file_resources input, and
        returns a DataFrame.

        Args:
            file_resources (dict): A dict with keys as filenames and values as
                full file path.
            blocksize (int):
        """
        raise NotImplementedError

    @classmethod
    def name(cls):
        return cls.__name__

    @staticmethod
    def list_databases():
        return DEFAULT_LIBRARIES

    def get_annotations(self, on: Union[str, List[str]],
                        columns: List[str],
                        agg: str = "unique",
                        keys: pd.Index = None):
        """Returns the Database's DataFrame such that it's indexed by :param
        index:, which then applies a groupby operation and aggregates all other
        columns by concatenating all unique values.

        Args:
            on (str): The column name of the DataFrame to join by.
            columns (list): a list of column names.
            agg (str): Function to aggregate when there is more than one values
                for each index instance. E.g. ['first', 'last', 'sum', 'mean',
                'size', 'concat'], default 'concat'.
            keys (pd.Index): The values on the `index` column to
                filter before performing the groupby-agg operations.

        Returns:
            DataFrame: A dataframe to be used for annotation
        """
        if not set(columns).issubset(set(self.data.columns)):
            raise Exception(
                f"The columns argument must be a list such that it's subset of the following columns in the dataframe. "
                f"These columns doesn't exist in database: {list(set(columns) - set(self.data.columns.tolist()))}"
            )

        # Select df columns including df. However, the `columns` list shouldn't contain the index column
        if on in columns:
            columns.pop(columns.index(on))

        select_col = columns + ([on] if not isinstance(on, list) else on)
        if isinstance(self.data, pd.DataFrame):
            df = self.data.filter(select_col, axis="columns")
        elif isinstance(self.data, dd.DataFrame):
            df = self.data[select_col]

        if keys is not None:
            if on == df.index.name:
                df = df.loc[keys]
            else:
                df = df[df[on].isin(keys)]

        # Groupby index
        if on != df.index.name and df.index.name in columns:
            groupby = df.reset_index().groupby(on)
        else:
            groupby = df.groupby(on)

        #  Aggregate by all columns by concatenating unique values
        if agg == "unique":
            if isinstance(df, pd.DataFrame):
                grouby_agg = groupby.agg({col: concat_uniques for col in columns})

            elif isinstance(df, dd.DataFrame):
                agg_func = dd.Aggregation(
                    name=agg,
                    chunk=lambda s1: s1.apply(concat_uniques),
                    agg=lambda s2: s2.obj,
                )
                grouby_agg = groupby.agg({col: agg_func for col in columns})

            else:
                raise Exception("Unsupported dataframe: {}".format(df))

        # Concatenate values into list
        elif agg == "concat":
            grouby_agg = groupby.agg({col: concat for col in columns})

        # Any other aggregation functions
        else:
            grouby_agg = groupby.agg({col: agg for col in columns})

        return grouby_agg

    def get_expressions(self, index):
        """
        Args:
            index:
        """
        return self.data.groupby(index).median(
        )  # TODO if index by gene, aggregate medians of transcript-level expressions



class Annotatable(ABC):
    """This abstract class provides an interface for the -omics
    (:class:`Expression`) to annotate its genes list with the external data
    downloaded from various databases. The database will be imported as
    attributes information to the genes's annotations, or interactions between
    the genes.
    """
    SEQUENCE_COL = "sequence"
    DISEASE_ASSOCIATIONS_COL = "disease_associations"

    def get_annotations(self):
        if hasattr(self, "annotations"):
            return self.annotations
        else:
            raise Exception(
                "{} must run initialize_annotations() first.".format(
                    self.name()))

    def get_annotation_expressions(self):
        if hasattr(self, "annotation_expressions"):
            return self.annotation_expressions
        else:
            raise Exception("{} must run annotate_expressions() first.".format(
                self.name()))

    def initialize_annotations(self, index=None):
        """
        Args:
            index:
            gene_list:
        """
        if index is None:
            index = self.get_genes_list()

        self.annotations: pd.DataFrame = pd.DataFrame(index=index)

    def annotate_attributes(self, database: Union[Database, pd.DataFrame], on: Union[str, List[str]],
                            columns: List[str], agg: str = "unique",
                            fuzzy_match: bool = False):
        """Performs a left outer join between the annotation and Database's
        DataFrame, on the index key. The index argument must be column present
        in both DataFrames. If there exists overlapping columns from the join,
        then .fillna() is used to fill NaN values in the old column with non-NaN
        values from the new column.

        Args:
            database (Database): Database which contains an dataframe.
            on (str): The column name which exists in both the annotations and
                Database dataframe to perform the join on.
            columns ([str]): a list of column name to join to the annotation.
            agg (str): Function to aggregate when there is more than one values
                for each index instance. E.g. ['first', 'last', 'sum', 'mean',
                'concat'], default 'concat'.
            fuzzy_match (bool): default False. Whether to join the annotation by applying a fuzzy match on the index
                with difflib.get_close_matches(). It can be slow and thus should only be used sparingly.
        """
        if not hasattr(self, "annotations"):
            raise Exception("Must run .initialize_annotations() on, ", self.__class__.__name__, " first.")

        if_pandas_df = isinstance(self.annotations, pd.DataFrame)
        old_index = self.annotations.index.name

        if isinstance(on, str) and on in self.annotations.columns:
            keys = self.annotations[on]
        elif isinstance(on, list) and set(on).issubset(self.annotations.columns):
            keys = self.annotations[on]
        elif on == self.annotations.index.name:
            keys = self.annotations.index
        else:
            keys = None

        if isinstance(database, pd.DataFrame):
            database_df = database.groupby(on).agg({col: concat_uniques for col in columns})
        else:
            database_df = database.get_annotations(on, columns=columns, agg=agg, keys=keys)

        if len(database_df) == 0:
            logging.warning("Database annotations is empty and has nothing to annotate.")
            return

        if fuzzy_match:
            database_df.index = database_df.index.map(
                lambda x: difflib.get_close_matches(x, self.annotations.index, n=1)[0])

        # performing join on the index column
        if isinstance(database_df, pd.DataFrame):
            new_annotations = self.annotations.join(database_df, on=on, rsuffix="_")
        # performing join on a different column
        elif isinstance(database_df, dd.DataFrame):
            new_annotations = dd.merge(self.annotations, database_df, how="left", on=on, suffixes=("_", ""))
            # Revert back to Pandas df if not previously a dask df
            if if_pandas_df and isinstance(new_annotations, dd.DataFrame):
                new_annotations = new_annotations.compute()

        # Merge columns if the database DataFrame has overlapping columns with existing column
        duplicate_cols = [col for col in new_annotations.columns if col[-1] == "_"]

        # Fill in null values then drop duplicate columns
        for new_col in duplicate_cols:
            old_col = new_col.strip("_")
            new_annotations[old_col] = new_annotations[old_col].fillna(new_annotations[new_col], axis=0)
            new_annotations = new_annotations.drop(columns=new_col)

        if old_index != self.annotations.index.name and old_index != self.annotations.index.names:
            self.annotations = self.annotations.set_index(old_index)

        # Assign the new results
        self.annotations = new_annotations

    def annotate_sequences(self,
                           database,
                           on: Union[str, List[str]],
                           agg="longest",
                           omic=None,
                           **kwargs):
        """Annotate a genes list (based on index) with a dictionary of
        <gene_name: sequence>. If multiple sequences per gene name, then perform
        some aggregation.

        Args:
            database (SequenceDatabase): The database
            on (str): The gene index column name.
            agg (str): The aggregation method, one of ["longest", "shortest", or
                "all"]. Default longest.
            omic (str): Default None. Declare the omic type to fetch sequences
                for.
            **kwargs:
        """
        if omic is None:
            omic = self.name()

        sequences = database.get_sequences(index=on, omic=omic, agg=agg, **kwargs)

        # Map sequences to the keys of `on` columns.
        if type(self.annotations.index) == pd.MultiIndex and self.annotations.index.names in on:
            seqs = self.annotations.index.get_level_values(on).map(sequences)
        elif self.annotations.index.name == on:
            seqs = self.annotations.index.map(sequences)
        elif isinstance(on, list):
            # Index is a multi columns
            seqs = pd.MultiIndex.from_frame(self.annotations.reset_index()[on]).map(sequences)
        else:
            seqs = pd.Index(self.annotations.reset_index()[on]).map(sequences)

        self.annotations[Annotatable.SEQUENCE_COL] = seqs

    def annotate_expressions(self, database, index, fuzzy_match=False):
        """Annotate :param database: :param index: :param fuzzy_match:

        Args:
            database:
            index:
            fuzzy_match:
        """
        self.annotation_expressions = pd.DataFrame(index=self.annotations.index)

        if self.annotations.index.name == index:
            self.annotation_expressions = self.annotation_expressions.join(
                database.get_expressions(index=index))
        else:
            raise Exception(f"index argument must be one of {database.data.index}")

    def annotate_interactions(self, database, index):
        """
        Args:
            database (Interactions):
            index (str):
        """
        raise NotImplementedError

    def annotate_diseases(self, database, on):
        """
        Args:
            database (DiseaseAssociation):
            on (str):
        """
        keys = self.annotations.index if on == self.annotations.index.name or on == self.annotations.index.names else \
        self.annotations[on]
        keys = pd.MultiIndex.from_frame(keys) if isinstance(keys, pd.DataFrame) else pd.Index(keys)
        disease_assocs = keys.map(database.get_disease_assocs(index=on, ))

        self.annotations[Annotatable.DISEASE_ASSOCIATIONS_COL] = disease_assocs

    def set_index(self, new_index):
        """Resets :param new_index: :type new_index: str

        Args:
            new_index:
        """
        self.annotations[new_index].fillna(self.annotations.index.to_series(),
                                           axis=0,
                                           inplace=True)
        self.annotations = self.annotations.reset_index().set_index(new_index)

    def get_rename_dict(self, from_index, to_index):
        """
        Utility function used to retrieve a lookup dictionary to convert from one index to
        another, e.g., gene_id to gene_name, obtained from two columns in the dataframe.

        Returns
            Dict[str, str]: the lookup dictionary.

        Args:
            from_index (str): an index on the DataFrame for key
            to_index:
        """
        if self.annotations.index.name in [from_index, to_index]:
            dataframe = self.annotations.reset_index()
        else:
            dataframe = self.annotations

        dataframe = dataframe[dataframe[to_index].notnull()]
        return pd.Series(dataframe[to_index].values,
                         index=dataframe[from_index]).to_dict()


DEFAULT_LIBRARIES = [
    "10KImmunomes"
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
    "TargetScan"
]


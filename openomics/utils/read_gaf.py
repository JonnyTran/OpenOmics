import gzip
from io import TextIOWrapper
from os.path import exists
from typing import List, Optional, Union, Dict, Callable

import dask.dataframe as dd
import pandas as pd
from Bio.UniProt.GOA import GAF20FIELDS, GAF10FIELDS, _gaf20iterator, _gaf10iterator
from filetype import filetype
from six import string_types
from six.moves import intern


def read_gaf(filepath_or_buffer, blocksize=None, compression: Optional[str] = None,
             column_converters: Dict[str, Callable] = None, usecols: List[str] = None, chunksize=1024 * 1024) -> Union[
    pd.DataFrame, dd.DataFrame]:
    """Parse a GTF into a dictionary mapping column names to sequences of
    values.

    Args:
        filepath_or_buffer (str or buffer object): Path to GTF file (may be gzip
            compressed) or buffer object such as StringIO
        blocksize (int): Number of blocksize for the dask dataframe. Default None to use pandas.DataFrame instead.
        compression (str): Compression type to be passed into dask.dataframe.read_table(). Default None.
        column_converters (dict, optional): Dictionary mapping column names to
            conversion functions. Will replace empty strings with None and
            otherwise passes them to given conversion function.
        usecols (list of str or None): Restrict which columns are loaded to the
            give set. If None, then load all columns.
    """
    if isinstance(filepath_or_buffer, string_types) and not exists(filepath_or_buffer):
        raise ValueError(f"GAF file does not exist: {filepath_or_buffer}")

    COLUMN_NAMES = infer_gaf_columns(filepath_or_buffer)

    if blocksize:
        assert isinstance(filepath_or_buffer, str), f'dd.read_table() must have `filepath_or_buffer` as a path, and ' \
                                                    f'if compressed, use the `compression` arg. ' \
                                                    f'`filepath_or_buffer`={filepath_or_buffer}'
        result_df = parse_gaf(filepath_or_buffer, column_names=COLUMN_NAMES,
                              blocksize=blocksize, compression=compression, chunksize=chunksize)
    else:
        result_df = parse_gaf(filepath_or_buffer, column_names=COLUMN_NAMES)

    if column_converters:
        for column_name, column_type in column_converters.items():
            result_df[column_name] = result_df[column_name].map(
                lambda s: column_type(s) if isinstance(s, str) and len(s) else None)

    if usecols is not None:
        column_names = set(result_df.columns)
        valid_columns = [c for c in usecols if c in column_names]
        result_df = result_df[valid_columns]

    return result_df


def parse_gaf(filepath_or_buffer, column_names=None, usecols=None,
              intern_columns=['DB', 'Evidence', 'Aspect', 'DB_Object_Type', 'Assigned_By'],
              list_dtype_columns=['Qualifier', 'DB:Reference', 'With', 'Synonym'],
              blocksize=None, compression=None, chunksize=1024 * 1024) \
    -> Union[pd.DataFrame, dd.DataFrame]:
    """

    Args:
        filepath_or_buffer ():
        column_names ():
        usecols ():
        intern_columns ():
        list_dtype_columns ():
        blocksize ():
        compression ():
        chunksize ():

    Returns:

    """

    def split_str(input: str, sep='|') -> List[str]:
        return input.split(sep)

    def parse_taxon(input: str) -> List[str]:
        taxons = input.split("|")
        return [s.strip("taxon:") for s in taxons]

    args = dict(
        sep="\t",
        comment="!",
        names=column_names,
        usecols=usecols,
        skipinitialspace=True,
        skip_blank_lines=True,
        on_bad_lines='error',
        engine="c",
        dtype={
            "Date": 'str',
            "Annotation_Extension": "str",
            "Gene_Product_Form_ID": "str",
        },
        date_parser=pd.to_datetime,
        converters={
            'Taxon_ID': parse_taxon,
            **{col: split_str for col in (list_dtype_columns if list_dtype_columns else [])}
        },
    )

    if blocksize:
        df = dd.read_table(filepath_or_buffer, compression=compression, blocksize=blocksize, **args)
        df['Date'] = dd.to_datetime(df['Date'])

    else:
        chunk_iterator = pd.read_table(filepath_or_buffer, chunksize=chunksize, parse_dates=['Date'], **args)
        dataframes = []
        try:
            for df in chunk_iterator:
                for intern_column in intern_columns:
                    df[intern_column] = df[intern_column].map(lambda s: intern(str(s)))

                dataframes.append(df)

            df = pd.concat(dataframes)
        except Exception as e:
            raise Exception("ParsingError:" + str(e))

    return df


def infer_gaf_columns(handle: Union[str, TextIOWrapper]) -> List[str]:
    """
    Grab first line of the file, filestream, or a compressed file, then determine the gaf version and return
    corresponding column names to expect from the tab-delimited .gaf file.
    Args:
        handle ():

    Returns:

    """
    if isinstance(handle, str) and exists(handle):
        filepath_ext = filetype.guess(handle)

        if filepath_ext is not None and filepath_ext.extension == 'gz':
            with gzip.open(handle, 'rt') as file:
                inline = file.readline()
        else:
            with open(handle, 'rt') as file:
                inline = file.readline()

    elif hasattr(handle, 'readline'):
        inline = handle.readline()
    else:
        raise IOError(f"`handle` of type {type(handle)} not supported.")

    if inline.strip().startswith("!gaf-version: 2"):
        return GAF20FIELDS
    elif inline.strip().startswith("!gaf-version: 1.0"):
        return GAF10FIELDS
    else:
        raise Exception(f"{inline} not supported.")


@DeprecationWarning
def gafiterator(handle: TextIOWrapper):
    inline = handle.readline()
    if inline.strip().startswith("!gaf-version: 2"):
        return _gaf20iterator(handle)
    elif inline.strip() == "!gaf-version: 1.0":
        return _gaf10iterator(handle)
    else:
        return _gaf20iterator(handle)

import gzip
import os.path
import warnings
from collections.abc import Iterable
from io import TextIOWrapper
from os.path import exists
from typing import List, Optional, Union, Dict, Callable

import dask.dataframe as dd
import numpy as np
import pandas as pd
from Bio.UniProt.GOA import GAF20FIELDS, GAF10FIELDS
from filetype import filetype
from logzero import logger
from six import string_types
from six.moves import intern


def read_gaf(filepath_or_buffer, index_col=None, keys=None, compression: Optional[str] = None,
             column_converters: Dict[str, Callable] = None, usecols: List[str] = None, chunksize=1024 * 1024,
             blocksize=None, ) \
    -> Union[pd.DataFrame, dd.DataFrame]:
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
    if isinstance(filepath_or_buffer, str) and "~" in filepath_or_buffer:
        filepath_or_buffer = os.path.expanduser(filepath_or_buffer)

    if isinstance(filepath_or_buffer, string_types) and not exists(filepath_or_buffer):
        raise ValueError(f"GAF file does not exist: {filepath_or_buffer}")

    COLUMN_NAMES = infer_gaf_columns(filepath_or_buffer)

    if blocksize:
        assert isinstance(filepath_or_buffer, str), f'dd.read_table() must have `filepath_or_buffer` as a path, and ' \
                                                    f'if compressed, use the `compression` arg. ' \
                                                    f'`filepath_or_buffer`={filepath_or_buffer}'
        result_df = parse_gaf(filepath_or_buffer, column_names=COLUMN_NAMES, index_col=index_col, keys=keys,
                              blocksize=blocksize,
                              compression=compression, chunksize=chunksize)
    else:
        result_df = parse_gaf(filepath_or_buffer, column_names=COLUMN_NAMES, index_col=index_col, keys=keys)

    if column_converters:
        for column_name, column_type in column_converters.items():
            result_df[column_name] = result_df[column_name].map(
                lambda s: column_type(s) if isinstance(s, str) and len(s) else None)

    if usecols is not None:
        column_names = set(result_df.columns)
        valid_columns = [c for c in usecols if c in column_names]
        result_df = result_df[valid_columns]

    return result_df


def parse_gaf(filepath_or_buffer, column_names=None, index_col=None, keys=None, usecols=None,
              intern_columns=['DB', 'Evidence', 'Aspect', 'DB_Object_Type', 'Assigned_By'],
              list_dtype_columns=['DB:Reference', 'With', 'Synonym'],
              blocksize=None, compression=None,
              chunksize=1024 * 1024, ) \
    -> Union[pd.DataFrame, dd.DataFrame]:
    """

    Args:
        filepath_or_buffer (str):
        column_names ():
        usecols ():
        intern_columns (): Only used when blocksize == None, i.e. when using Pandas DataFrames
        list_dtype_columns ():
        blocksize ():
        compression ():
        chunksize ():

    Returns:

    """

    def split_str(input: str, sep='|') -> Optional[np.ndarray]:
        if isinstance(input, str):
            return np.array(input.split(sep))
        elif isinstance(input, Iterable):
            return input
        else:
            return None

    def parse_taxon(input: str, sep='|') -> Optional[np.ndarray]:
        if isinstance(input, str):
            return np.array([s.replace("taxon:", '').strip() for s in input.split(sep)])
        elif isinstance(input, Iterable):
            return input
        else:
            return None

    parse_args = dict(
        sep="\t",
        comment="!",
        names=column_names,
        usecols=usecols,
        skipinitialspace=True,
        skip_blank_lines=True,
        on_bad_lines='error',
        dtype='str',
    )

    if blocksize:
        if filepath_or_buffer.endswith('.parquet'):
            # Raw GAF file in .parquet format
            df: dd.DataFrame = dd.read_parquet(filepath_or_buffer, usecols=usecols, chunksize=blocksize)
        else:
            df: dd.DataFrame = dd.read_table(filepath_or_buffer, compression=compression,
                                             blocksize=blocksize if blocksize > 10 else None, **parse_args)

        if index_col:
            # Filter
            if keys is not None and df.index.name != index_col:
                df = df.loc[df[index_col].isin(keys)]
            elif keys is not None and df.index.name == index_col:
                df = df.loc[df.index.isin(keys)]

            # Set index
            if df.index.name != index_col:
                logger.info(f"Setting index at {index_col}, existing index {df.index.name}")
                df = df.set_index(index_col, sorted=False).persist()

            # Compute division of partitions in index
            if not df.known_divisions:
                df.divisions = df.compute_current_divisions()

        if not 'processed' in filepath_or_buffer:
            # Transform columns
            assign_fn = {col: df[col].map(split_str) for col in df.columns.intersection(list_dtype_columns)}
            if 'Taxon_ID' in df.columns:
                assign_fn['Taxon_ID'] = df['Taxon_ID'].map(parse_taxon)
            if 'Date' in df.columns and df['Date'].dtype == 'O':
                assign_fn['Date'] = dd.to_datetime(df['Date'])
            df = df.assign(**assign_fn)

    else:
        parse_args['converters'] = {
            'Taxon_ID': parse_taxon,
            **{col: split_str for col in (list_dtype_columns if list_dtype_columns else [])}
        }
        chunk_iterator = pd.read_table(filepath_or_buffer, chunksize=chunksize, index_col=index_col,
                                       parse_dates=['Date'], date_parser=pd.to_datetime, **parse_args)
        dataframes = []
        try:
            for df in chunk_iterator:
                for intern_column in intern_columns:
                    df[intern_column] = df[intern_column].map(lambda s: intern(str(s)))

                if keys is not None and df.index.name != index_col:
                    df = df.loc[df[index_col].isin(keys)]
                elif keys is not None and df.index.name == index_col:
                    df = df.loc[df.index.isin(keys)]

                dataframes.append(df)

            df = pd.concat(dataframes)
        except Exception as e:
            raise Exception("ParsingError:" + str(e))

    return df


def infer_gaf_columns(filepath_or_buffer: Union[str, TextIOWrapper], default_gaf_fields=GAF20FIELDS) -> List[str]:
    """
    Grab first line of the file, filestream, or a compressed file, then determine the gaf version and return
    corresponding column names to expect from the tab-delimited .gaf file.
    Args:
        filepath_or_buffer ():

    Returns:

    """
    if isinstance(filepath_or_buffer, str) and exists(filepath_or_buffer) and os.path.isfile(filepath_or_buffer):
        filepath_ext = filetype.guess(filepath_or_buffer)

        if filepath_ext is not None and filepath_ext.extension == 'gz':
            with gzip.open(filepath_or_buffer, 'rt') as file:
                inline = file.readline()
        else:
            with open(filepath_or_buffer, 'rt') as file:
                inline = file.readline()

    elif hasattr(filepath_or_buffer, 'readline'):
        inline = filepath_or_buffer.readline()
    else:
        warnings.warn(f"`filepath_or_buffer`={filepath_or_buffer} not supported to peek first line and "
                      f"infer GAF version for column names. Defaulting to `GAF20FIELDS`")
        return default_gaf_fields

    if inline.strip().startswith("!gaf-version: 2"):
        return GAF20FIELDS
    elif inline.strip().startswith("!gaf-version: 1.0"):
        return GAF10FIELDS
    else:
        raise Exception(f"{inline} not supported.")



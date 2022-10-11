from collections.abc import Iterable
from typing import Union, Dict, List, Tuple, Optional, Set

import numpy as np
import pandas as pd
from dask import dataframe as dd
from logzero import logger


def has_iterables(series: Union[dd.Series, pd.Series], n=10) -> bool:
    """
    Check whether any element in series is an Iterable, e.g. list, set, or np.ndarray.
    Args:
        series ():
        n (int): number of elements to test

    Returns:
        bool
    """
    if isinstance(series, (dd.Series, dd.Index)):
        values = series.head(n=n, npartitions=-1)
    elif isinstance(series, pd.Series):
        values = series.head(n=n)
    elif isinstance(series, pd.Index):
        values = series[:n]
    else:
        return False

    is_iterables = values.map(lambda x: not isinstance(x, str) and isinstance(x, Iterable)).any()
    return is_iterables


def match_iterable_keys(left: Union[dd.Series, pd.Series], right: Union[dd.Series, pd.Series]) \
    -> Tuple[dd.Series, dd.Series]:
    left_iterables = has_iterables(left)
    right_iterables = has_iterables(right)

    def _list_to_key(list_values: List[str], possibilities: Set[str]) -> Optional[str]:
        if list_values is None:
            return None

        for key in list_values:
            if key in possibilities:
                return key

    def _list_to_list(list_values: List[str], multi_possibilities: List[Set[str]]) -> Optional[str]:
        if list_values is None:
            return None

        for possibilities in multi_possibilities:
            if isinstance(possibilities, set) and list_values is not None:
                match = possibilities.intersection(list_values)
                if len(match):
                    return "|".join(sorted(match))

    left_on, right_on = left, right

    if left_iterables and not right_iterables:
        possibilities = set(right.dropna()) if isinstance(right, (pd.Series, pd.Index)) else set(
            right.dropna().compute())

        left_on = left.map(lambda values: _list_to_key(values, possibilities))
    elif not left_iterables and right_iterables:
        possibilities = set(left.dropna()) if isinstance(left, (pd.Series, pd.Index)) else set(
            left.dropna().compute())

        right_on = right.map(lambda values: _list_to_key(values, possibilities))
    elif left_iterables and right_iterables:
        right_possibilities = right.map(lambda x: set(x) if isinstance(x, Iterable) and not isinstance(x, str) else x)
        if isinstance(right_possibilities, (dd.Series, dd.Index)):
            right_possibilities = right_possibilities.compute()
        left_on = left.map(lambda values: _list_to_list(values, right_possibilities))

        left_possibilities = left.map(lambda x: set(x) if isinstance(x, Iterable) and not isinstance(x, str) else x)
        if isinstance(left_possibilities, (dd.Series, dd.Index)):
            left_possibilities = left_possibilities.compute()
        right_on = right.map(lambda values: _list_to_list(values, left_possibilities))

    return left_on, right_on


def drop_duplicate_columns(df: Union[pd.DataFrame, dd.DataFrame]) -> Union[pd.DataFrame, dd.DataFrame]:
    """
    Args:
        df:
    """
    if df.columns.duplicated().any():
        _, i = np.unique(df.columns, return_index=True)
        df = df.iloc[:, i]

    return df


def filter_rows(df: pd.DataFrame, filters: Union[str, Dict[str, List]], case: bool = True):
    """

    Args:
        df (pd.DataFrame):
        filters (str or dict):
            Either a pandas query expression or a dict of column names for keys and matching values.
        case (bool): Default True.
            Whether to match case in pd.Series.str.contains if filters is a dict of values.

    Returns:

    """
    num_samples = df.shape[0]
    if filters is None:
        return df

    elif isinstance(filters, str):
        df = df.query(filters)

    elif isinstance(filters, dict):
        for col, values in filters.items():
            if col not in df.columns:
                logger.warn("Filter key `", col, "` must be in one of ", df.columns)
                continue

            if isinstance(values, list):
                values = {val.upper() for val in values}
                if case is True:
                    df = df.loc[df[col].isin(values)]
                else:
                    df = df.loc[df[col].str.upper().isin(values)]
            elif isinstance(values, str):
                df = df.loc[df[col].str.contains(values, case=case)]
            else:
                df = df.loc[df[col] == values]

    if isinstance(df, pd.DataFrame) and df.shape[0] == 0:
        logger.warn(f"Dataframe is empty ({df.shape}) because of query: {filters}")

    if isinstance(num_samples, int):
        logger.info(f'Removed {num_samples - df.shape[0]} rows from query: {filters}')

    return df

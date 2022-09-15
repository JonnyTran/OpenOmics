from collections import defaultdict
from collections.abc import Iterable
from typing import Union, List, Dict, Callable, Mapping

import dask.dataframe as dd
import numpy as np
import pandas as pd


def get_aggregator(keyword: str, use_dask=False) -> Union[str, Callable, dd.Aggregation]:
    """

    Args:
        keyword ():
        use_dask ():

    Returns:

    """
    if keyword == "unique":
        if use_dask:
            func = dd.Aggregation(
                name=keyword,
                chunk=lambda s1: s1.apply(concat_uniques),
                agg=lambda s2: s2.obj,
            )
        else:
            func = concat_uniques

    elif keyword == "concat":
        # Concatenate values into list
        func = concat
    else:
        # Any other aggregation keywords or callable function
        func = keyword

    return func


def get_multi_aggregators(agg: str, agg_for: Dict[str, Union[str, Callable, dd.Aggregation]] = None, use_dask=False) \
    -> Mapping[str, Union[str, dd.Aggregation]]:
    """

    Args:
        agg ():
        agg_for ():
        use_dask ():

    Returns:

    """
    if agg_for is None:
        agg_for = {}

    agg_for = {col: get_aggregator(keyword, use_dask=use_dask) for col, keyword in agg_for.items()}
    col_aggregators = defaultdict(get_aggregator(agg, use_dask=use_dask), agg_for)

    return col_aggregators


def concat_uniques(series: pd.Series) -> Union[str, List, np.ndarray, None]:
    """ An aggregation custom function to be applied to each column of a groupby
    Args:
        series (pd.Series): Entries can be either a string or a list of strings.
    Returns:
        unique_values
    """
    series = series.dropna()
    if series.empty:
        return None

    is_str_idx = series.map(type) == str

    if series.map(lambda x: isinstance(x, Iterable)).any():
        if (is_str_idx).any():
            # Convert mixed dtypes to lists
            series.loc[is_str_idx] = series.loc[is_str_idx].map(lambda s: [s] if len(s) else None)
        return np.unique(np.hstack(series))

    elif is_str_idx.any():
        concat_str = series.astype(str).unique()
        if len(concat_str):  # Avoid empty string
            return concat_str

    else:
        return series.tolist()


def concat(series: pd.Series) -> Union[str, List, np.ndarray, None]:
    """
    Args:
        series (pd.Series): Entries can be either a string or a list of strings.
    """
    series = series.dropna()
    if series.empty:
        return

    is_str_idx = series.map(type) == str
    if series.map(lambda x: isinstance(x, Iterable)).any():
        if (is_str_idx).any():
            # Convert mixed dtypes to lists
            series.loc[is_str_idx] = series.loc[is_str_idx].map(lambda s: [s] if len(s) else None)
        return np.hstack(series)

    elif is_str_idx.any():
        concat_str = series.astype(str).tolist()
        if len(concat_str):  # Avoid empty string
            return concat_str

    else:
        return series.tolist()

def drop_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Args:
        df:
    """
    if df.columns.duplicated().any():
        _, i = np.unique(df.columns, return_index=True)
        df = df.iloc[:, i]

    return df

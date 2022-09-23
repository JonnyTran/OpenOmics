from collections import defaultdict
from typing import Union, List, Dict, Callable, Mapping, Iterable

import dask.dataframe as dd
import numpy as np
import pandas as pd
from pandas.core.groupby import SeriesGroupBy


def get_agg_func(keyword: str, use_dask=False) -> Union[str, Callable, dd.Aggregation]:
    """

    Args:
        keyword (str):
        use_dask (bool): Whether to create a dd.Aggregation

    Returns:
        func (callable): a callable function, pandas aggregator func name, or a Dask Aggregation.
    """
    if keyword == "unique" and use_dask:
        # get unique values (in a list-like np.array) from each groupby key
        func = concat_unique_dask_agg()
    elif keyword == "unique" and not use_dask:
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

    col2func = {col: get_agg_func(keyword, use_dask=use_dask) for col, keyword in agg_for.items()}
    col_aggregators = defaultdict(lambda: get_agg_func(agg, use_dask=use_dask), col2func)

    return col_aggregators


def concat_unique_dask_agg() -> dd.Aggregation:
    def chunk(s: pd.Series) -> pd.Series:
        '''
        The function applied to the individual partition (map)
        '''

        def to_list(x: Union[str, List, np.ndarray]):
            if isinstance(x, str):
                return [x]
            elif isinstance(x, np.ndarray):
                return x
            elif isinstance(x, Iterable):
                if any(isinstance(a, Iterable) for a in x):
                    return list(set(np.hstack(x)))
                else:
                    return list(set(x))
            else:
                return [x]

        return s.apply(to_list)

    def agg(s: SeriesGroupBy) -> pd.Series:
        '''
        The function which will aggregate the result from all the partitions(reduce)
        '''
        s = s._selected_obj
        return s.groupby(level=list(range(s.index.nlevels))).apply(np.hstack)

    def finalize(s) -> pd.Series:
        '''
        The optional functional that will be applied to the result of the agg_tu functions
        '''
        return s.apply(lambda arr: np.unique(arr[~pd.isna(arr)]))

    func = dd.Aggregation(
        name='unique',
        chunk=chunk,
        agg=agg,
        finalize=finalize
    )
    return func


def merge_values(a: Union[str, None, Iterable], b: Union[str, None, Iterable]) -> Union[np.ndarray, str, None]:
    """
    Used as function in pd.combine() or dd.combine()
    Args:
        a (Union[str,None,Iterable]):
        b (Union[str,None,Iterable]):

    Returns:
        combined_value (Union[np.ndarray, str, None])
    """
    a_isna = pd.isna(a)
    b_isna = pd.isna(b)
    if a_isna is True or (isinstance(a_isna, Iterable) and all(a_isna)):
        return b
    elif b_isna is True or (isinstance(b_isna, Iterable) and all(b_isna)):
        return a
    elif isinstance(a, str) and isinstance(b, str):
        return np.array([a, b])
    elif not isinstance(a, Iterable) and isinstance(b, Iterable):
        return np.hstack([[a], b])
    elif isinstance(a, Iterable) and not isinstance(b, Iterable):
        return np.hstack([a, [b]])
    elif isinstance(a, Iterable) and isinstance(b, Iterable):
        return np.hstack([a, b])
    else:
        return b


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

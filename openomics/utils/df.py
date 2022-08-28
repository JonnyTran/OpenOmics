from typing import Union, List

import numpy as np
import pandas as pd


def concat_uniques(series: pd.Series, sep="|") -> Union[str, List, np.ndarray, None]:
    """ An aggregation custom function to be applied to each column of a groupby
    Args:
        series (pd.Series):
    """
    series = series.dropna()
    if series.empty:
        return None

    is_str_idx = series.map(type) == str

    if series.map(lambda x: isinstance(x, (list, tuple))).any():
        if (is_str_idx).any():
            # Convert mixed dtypes to lists
            series.loc[is_str_idx] = series.loc[is_str_idx].map(lambda s: [s] if len(s) else None)
        return np.unique(np.hstack(series))

    elif is_str_idx.any():
        concat_str = sep.join(series.astype(str).unique())
        if len(concat_str):
            return concat_str

    else:
        return series.tolist()


def concat(series: pd.Series, sep="|") -> Union[str, List, np.ndarray, None]:
    """
    Args:
        series (pd.Series):
    """
    series = series.dropna()
    if series.empty:
        return

    is_str_idx = series.map(type) == str
    if series.map(lambda x: isinstance(x, (list, tuple))).any():
        if (is_str_idx).any():
            # Convert mixed dtypes to lists
            series.loc[is_str_idx] = series.loc[is_str_idx].map(lambda s: [s] if len(s) else None)
        return np.hstack(series)

    elif is_str_idx.any():
        concat_str = sep.join(series.astype(str))
        if len(concat_str):
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

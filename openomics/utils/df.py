from typing import Union, List

import numpy as np
import pandas as pd


def concat_uniques(series: pd.Series) -> Union[str, List, np.ndarray]:
    """ An aggregation custom function to be applied to each column of a groupby
    Args:
        series (pd.Series):
    """
    if series.empty:
        return None
    series = series.dropna()

    if series.map(lambda x: isinstance(x, (list, tuple))).any():
        return np.unique(np.hstack(series))
    else:
        return "|".join(series.astype(str).unique())


def concat(series: pd.Series) -> Union[str, List, np.ndarray]:
    """
    Args:
        series (pd.Series):
    """
    if series.empty:
        return None
    series = series.dropna()

    if series.map(lambda x: isinstance(x, (list, tuple))).any():
        return np.hstack(series)
    else:
        return "|".join(series.astype(str))


def drop_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Args:
        df:
    """
    if df.columns.duplicated().any():
        _, i = np.unique(df.columns, return_index=True)
        df = df.iloc[:, i]

    return df

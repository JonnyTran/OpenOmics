from typing import Union, Dict, List

import numpy as np
import pandas as pd
from dask import dataframe as dd
from logzero import logger


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
        logger.info(f"Dataframe is empty ({df.shape}) because of query: {filters}")
    return df

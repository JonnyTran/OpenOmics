from typing import Union, Dict, List

import numpy as np
import pandas as pd
from dask import dataframe as dd


def drop_duplicate_columns(df: Union[pd.DataFrame, dd.DataFrame]) -> Union[pd.DataFrame, dd.DataFrame]:
    """
    Args:
        df:
    """
    if df.columns.duplicated().any():
        _, i = np.unique(df.columns, return_index=True)
        df = df.iloc[:, i]

    return df


def filter_rows(df: pd.DataFrame, filters: Dict[str, List], case: bool = False):
    if filters is None:
        return df

    for col, values in filters.items():
        if col not in df.columns:
            print("Filter key `", col, "` must be in one of ", df.columns)
            continue
        n_rows = df.shape[0]

        if isinstance(values, list):
            if case is False:
                df = df.loc[df[col].str.upper().isin([val.upper() for val in values])]
            else:
                df = df.loc[df[col].isin(values)]
        elif isinstance(values, str):
            df = df.loc[df[col].str.contains(values, case=case)]
        else:
            df = df.loc[df[col] == values]

        print("INFO: Removed ", n_rows - df.shape[0], " rows with `", col, "` != ", values)

    assert df.shape[0] > 0, f"ERROR: Dataframe is empty ({df.shape}) because of filter: {filters}"
    return df

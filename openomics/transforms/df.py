from typing import Union

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

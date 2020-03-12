import pandas as pd


def concat_uniques(series: pd.Series):
    series = series.dropna().astype(str)
    if not series.empty:
        return "|".join(series.unique())
    else:
        return None


def concat(series: pd.Series):
    series = series.dropna().astype(str)
    if not series.empty:
        return "|".join(series)
    else:
        return None

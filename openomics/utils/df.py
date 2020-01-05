

def concat_uniques(x):
    agg = x.dropna().astype(str)
    if not agg.empty:
        return "|".join(agg.unique())
    else:
        return None

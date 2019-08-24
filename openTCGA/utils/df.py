

def concat_uniques_agg(x):
    agg = x.dropna()
    if not agg.empty:
        return "|".join(agg.unique())
    else:
        return None
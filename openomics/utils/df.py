
def concat_uniques(df_column):
    before_shape = df_column.shape
    df_column = df_column.dropna().astype(str)
    assert before_shape == df_column.shape

    if not df_column.empty:
        return "|".join(df_column.unique())
    else:
        return None

def concat(df_column):
    df_column = df_column.dropna().astype(str)
    if not df_column.empty:
        return "|".join(df_column)
    else:
        return None

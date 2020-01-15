

def concat_uniques(df_column):
    df_column = df_column.dropna().astype(str)
    if not df_column.empty:
        return "|".join(df_column.unique())
    else:
        return None

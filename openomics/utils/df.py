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


def slice_adj(adj, node_list: list, nodes_A, nodes_B=None):
    if nodes_B is None:
        idx = [node_list.index(node) for node in nodes_A]
        return adj[idx, :][:, idx]
    else:
        idx_A = [node_list.index(node) for node in nodes_A]
        idx_B = [node_list.index(node) for node in nodes_B]
        return adj[idx_A, :][:, idx_B]

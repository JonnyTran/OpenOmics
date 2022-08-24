from typing import Union, List, Dict, Tuple

import networkx as nx


def to_scipy_adjacency(g: nx.DiGraph, nodes: Union[List[str], Dict[str, List[str]]],
                       edge_types: Union[List[str], Tuple[str, str, str]] = None,
                       reverse=False,
                       format="coo", d_ntype="_N"):
    if reverse:
        g = g.reverse(copy=True)

    if not isinstance(g, nx.MultiGraph):
        raise NotImplementedError

    if not isinstance(edge_types, (list, tuple, set)):
        edge_types = ["_E"]

    edge_index_dict = {}
    for etype in edge_types:
        if isinstance(g, nx.MultiGraph) and isinstance(etype, str):
            edge_subgraph = g.edge_subgraph([(u, v, e) for u, v, e in g.edges if e == etype])
            nodes_A = nodes
            nodes_B = nodes
            metapath = (d_ntype, etype, d_ntype)

        elif isinstance(g, nx.MultiGraph) and isinstance(etype, tuple) and isinstance(nodes, dict):
            metapath: Tuple[str, str, str] = etype
            head, etype, tail = metapath
            edge_subgraph = g.edge_subgraph([(u, v, e) for u, v, e in g.edges if e == etype])

            nodes_A = nodes[head]
            nodes_B = nodes[tail]

        elif etype == "_E":
            edge_subgraph = g.edges
            nodes_A = nodes
            nodes_B = nodes
            metapath = (d_ntype, etype, d_ntype)
        else:
            raise Exception(f"Edge types `{edge_types}` is ill formed.")

        biadj = nx.bipartite.biadjacency_matrix(edge_subgraph, row_order=nodes_A, column_order=nodes_B,
                                                format="coo")

        if format == "coo":
            edge_index_dict[metapath] = (biadj.row, biadj.col)
        elif format == "pyg":
            import torch
            edge_index_dict[metapath] = torch.stack(
                [torch.tensor(biadj.row, dtype=torch.long),
                 torch.tensor(biadj.col, dtype=torch.long)])
        else:
            edge_index_dict[metapath] = biadj

    return edge_index_dict


def slice_adj(adj, node_list: list, nodes_A, nodes_B=None):
    """
    Args:
        adj:
        node_list (list):
        nodes_A:
        nodes_B:
    """
    if nodes_B is None:
        idx = [node_list.index(node) for node in nodes_A]
        return adj[idx, :][:, idx]
    else:
        idx_A = [node_list.index(node) for node in nodes_A]
        idx_B = [node_list.index(node) for node in nodes_B]
        return adj[idx_A, :][:, idx_B]

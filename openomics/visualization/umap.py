import matplotlib.pyplot as plt
import numpy as np
import umap
from mpl_toolkits.mplot3d import Axes3D


def d3_umap(X, y_km, heat=None):
    fig = plt.figure(figsize=(6, 6))
    ax = Axes3D(fig)

    reducer = umap.UMAP(random_state=1234, n_components=3)
    X_embedded = reducer.fit_transform(X)

    node_colors = get_node_colormap(y_km)

    x, y, z = X_embedded[:, 0], X_embedded[:, 1], X_embedded[:, 2]
    ax.scatter(x, y, z, marker='o', c=node_colors, cmap=plt.cm.get_cmap('Set1', 5))

    ax.set_xlabel('UMAP-X')
    ax.set_ylabel('UMAP-Y')
    ax.set_zlabel('UMAP-Z')
    plt.show()
    return reducer


def get_node_colormap(node_label):
    if type(node_label) == list:
        node_labels = node_label
        sorted_node_labels = sorted(set(node_labels), reverse=True)
        colors = np.linspace(0, 1, len(sorted_node_labels))
        node_colormap = {f: colors[sorted_node_labels.index(f)] for f in set(node_labels)}
        node_colors = [node_colormap[n] if n in node_colormap.keys() else None for n in node_labels]

    elif node_label.dtype == "object":
        node_labels = node_label.str.split("|", expand=True)[0]
        sorted_node_labels = sorted(node_labels.unique(), reverse=True)
        colors = np.linspace(0, 1, len(sorted_node_labels))
        node_colormap = {f: colors[sorted_node_labels.index(f)] for f in node_labels.unique()}
        node_colors = [node_colormap[n] if n in node_colormap.keys() else None for n in node_labels]

    elif node_label.dtype == "float":
        node_labels = node_label.values
        node_colormap = None
        node_colors = [n / node_labels.max() for n in node_labels]
    return node_colors

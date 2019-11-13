import matplotlib.pyplot as plt
import umap
from mpl_toolkits.mplot3d import Axes3D


def d3_umap(X, y_km, heat=None):
    fig = plt.figure(figsize=(6, 6))
    ax = Axes3D(fig)

    reducer = umap.UMAP(random_state=1234, n_components=3)
    X_embedded = reducer.fit_transform(X)
    x, y, z = X_embedded[:, 0], X_embedded[:, 1], X_embedded[:, 2]
    ax.scatter(x, y, z, marker='o', c=y_km, cmap=plt.cm.get_cmap('Set1', 5))

    ax.set_xlabel('UMAP-X')
    ax.set_ylabel('UMAP-Y')
    ax.set_zlabel('UMAP-Z')
    plt.show()
    return reducer

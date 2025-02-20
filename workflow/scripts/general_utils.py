import sklearn
import numpy as np
import pandas as pd


def get_knn_labels(adata, label_key, obsm=None, knnidx=None, n_neighbors=None, radius=None, n_jobs=-1):
    """
    Compute the k-nearest neighbors (KNN) labels for each cell in the anndata object.

    Parameters
    ----------
    adata : AnnData
        The anndata object containing the data to be labeled.
    label_key : str
        The key in adata.obs containing the true labels.
    obsm : str
        The key in adata.obsm containing the KNN data.
    knnidx : array-like, optional
        The indices of the nearest neighbors. If not provided, it is computed from the knn_key.
    n_neighbors : int, optional
        The number of nearest neighbors to consider. If not provided, it is computed from the knn_key.
    radius : float, optional
        The radius of the ball to consider neighbors in. If not provided, it is computed from the knn_key.
    n_jobs : int, optional
        The number of jobs to run in parallel. If not provided, it is set to -1, which means all available cores.


    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the KNN labels for each cell in the anndata object.
    """
    if knnidx is None:
        nn = sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors, radius=radius, n_jobs=n_jobs).fit(
            adata.obsm[obsm]
        )
        if n_neighbors is not None and radius is not None:
            raise ValueError("Please provide either n_neighbors or radius, but not both.")
        if n_neighbors is not None:
            knndis, knnidx = nn.kneighbors(adata.obsm[obsm])
        elif radius is not None:
            knndis, knnidx = nn.radius_neighbors(adata.obsm[obsm])

    df_dummies = pd.get_dummies(adata.obs[label_key])

    # Convert df_dummies to a numpy array for efficient indexing
    if not isinstance(df_dummies, np.ndarray):
        dummy_array = np.array(df_dummies)

    if knnidx.ndim == 2:
        # knnidx contains same length list of neighbor indices for each donor
        knnlabels = dummy_array[knnidx].sum(1)
    else:
        # Initialize an empty list to store the summed labels
        knnlabels = []

        # Loop over each row in knnidx
        for neighbors in knnidx:
            # Get the one-hot encoded labels for the current neighbors
            neighbor_labels = dummy_array[neighbors]

            # Sum the labels across the neighbors (axis=0 sums column-wise)
            summed_labels = neighbor_labels.sum(axis=0)

            # Append the summed labels to the list
            knnlabels.append(summed_labels)

        # Convert the list back to a numpy array (optional)
        knnlabels = np.array(knnlabels)

    knnlabels = pd.DataFrame(knnlabels, index=adata.obs.index, columns=df_dummies.columns)
    return knnlabels

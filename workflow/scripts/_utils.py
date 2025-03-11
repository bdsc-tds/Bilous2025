import sklearn
import numpy as np
import pandas as pd
import scipy
import gseapy
from sklearn.model_selection import train_test_split, permutation_test_score, cross_validate
from sklearn.linear_model import LogisticRegression
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import StandardScaler


def get_knn_labels(
    adata, label_key, obsm=None, knnidx=None, n_neighbors=None, radius=None, n_jobs=-1, return_sparse_neighbors=False
):
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
    return_sparse_neighbors : bool, optional
        Return neighbors matrix as sparse in addition to knndis, knnidx list

    Returns
    -------
    knnlabels: pd.DataFrame
        A pandas DataFrame containing the KNN labels for each cell in the anndata object.
    knndis, knnidx: np.array
        Distance and indices of nearest neighbors
    knn_graph:
        if return_neighbors_sparse, additional return of sparse kNN matrix
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

    if return_sparse_neighbors:
        knn_graph = nn.kneighbors_graph()
        return knnlabels, knndis, knnidx, knn_graph
    else:
        return knnlabels, knndis, knnidx


def logit_pvalue(model, x):
    """Calculate z-scores for scikit-learn LogisticRegression.
    parameters:
        model: fitted sklearn.linear_model.LogisticRegression with intercept and large C
        x:     matrix on which the model was fit
    This function uses asymtptics for maximum likelihood estimates.
    """
    p = model.predict_proba(x)
    n = len(p)
    m = len(model.coef_[0]) + 1
    coefs = np.concatenate([model.intercept_, model.coef_[0]])
    x_full = np.matrix(np.insert(np.array(x), 0, 1, axis=1))
    ans = np.zeros((m, m))
    for i in range(n):
        ans = ans + np.dot(np.transpose(x_full[i, :]), x_full[i, :]) * p[i, 1] * p[i, 0]
    vcov = np.linalg.inv(np.matrix(ans))
    se = np.sqrt(np.diag(vcov))
    t = coefs / se
    p = (1 - scipy.stats.norm.cdf(abs(t))) * 2
    return p


def get_marker_rank_significance(rnk, gene_set, top_n=None):
    """
    Calculate the significance of marker ranks using pre-ranked GSEA and hypergeometric test.
    rnk must contain all features, its length is assumed to be the total number N for hypergoemetric testing

    Parameters:
    - rnk (pd.Series): index with feature names and values representing feature scores with decreasing ranks.
    - adata (AnnData): Anndata object containing the gene expression data.
    - gene_set (list): List of genes to test for enrichment.
    - top_n (optional, int): Number of top-ranked genes to consider for the hypergeometric test.

    Returns:
    - markers_rank_significance (pd.DataFrame): DataFrame containing the significance results.
    """

    if np.all(rnk == 0.0):
        print("Warning: all zero importance vector. Returning empty dataframe.")
        return pd.DataFrame(
            np.nan,
            index=[0],
            columns=[
                "Name",
                "Term",
                "ES",
                "NES",
                "NOM p-val",
                "FDR q-val",
                "FWER p-val",
                "Tag %",
                "Gene %",
                "Lead_genes",
                "hypergeometric_pvalue",
            ],
        )
    # Calculate marker rank significance from pre-ranked GSEA
    markers_rank_significance = gseapy.prerank(
        rnk=rnk,
        gene_sets=[{"gene_set": gene_set}],
        min_size=0,
    ).res2d

    # Calculate marker rank significance from hypergeometric test
    if top_n is not None:
        N = len(rnk)  # Total genes in ranked list
        K = len(gene_set)  # Genes in the pathway/set of interest
        x = np.isin(rnk.index[:top_n], gene_set).sum()  # Overlapping genes in top n

        # Add hypergeometric p-value to the results
        markers_rank_significance["hypergeometric_pvalue"] = scipy.stats.hypergeom.sf(x - 1, N, K, top_n)
        markers_rank_significance[f"n_hits_{top_n=}"] = x

    return markers_rank_significance


def logreg(
    X,
    y,
    feature_names=None,
    scoring="f1",
    test_size=0.2,
    n_permutations=30,
    n_repeats=5,
    max_iter=1000,
    random_state=0,
    importance_mode="coef",
):
    """
    Perform logistic regression with permutation test and compute feature importances.

    Parameters:
    - X (array-like): Input data for model training/testing.
    - y (vector-like): Input vector of labels for model training/testing.
    - feature_names (vector-like): Names of X features (optional).
    - scoring (str): Scoring metric for the permutation test (e.g., 'f1', 'accuracy').
    - test_size (float): Proportion of data to use for testing.
    - max_iter (int): Maximum number of iterations for the logistic regression model.
    - n_permutations (int): Number of permutations for the permutation test.
    - n_repeats (int): Number of repeats for the permutation importance calculation.
    - random_state (int): Random seed for reproducibility.
    - importance_mode (str): Mode for feature importance calculation ('permutation' or 'coef').

    Returns:
    - df_permutations (pd.DataFrame): Summary of permutation test results.
    - df_importances (pd.DataFrame): Feature importances from permutation importance.
    """

    # Initialize logistic regression model
    model = LogisticRegression(max_iter=max_iter)

    # Empirical p-value calculation using permutation test
    score, perm_scores, p_value = permutation_test_score(
        model, X, y, scoring=scoring, n_permutations=n_permutations, n_jobs=-1, verbose=1
    )

    # Summarize permutation test results
    df_permutations = pd.DataFrame(
        [[score, perm_scores.mean(), perm_scores.std(), p_value]],
        columns=[f"{scoring}_score", f"perm_mean{scoring}_score", f"perm_std{scoring}_score", "p_value"],
    )
    df_permutations["effect_size"] = (
        df_permutations[f"{scoring}_score"] - df_permutations[f"perm_mean{scoring}_score"]
    ) / df_permutations[f"perm_std{scoring}_score"]

    # Fit the model and compute feature importances from permutations

    if importance_mode == "permutation":
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, stratify=y, random_state=random_state
        )
        model.fit(X_train, y_train)
        importances = permutation_importance(
            model, pd.DataFrame.sparse.from_spmatrix(X_test), y_test, scoring=scoring, n_repeats=n_repeats, n_jobs=-1
        )
        importances.pop("importances")
        df_importances = pd.DataFrame(importances, index=feature_names).sort_values("importances_mean", ascending=False)

    elif importance_mode == "coef":
        model.fit(X, y)
        # Feature importances from model coefs
        # cv_results = cross_validate(model,X,y,return_estimator=True, scoring=scoring, n_jobs=-1)
        # importances = np.std(X, axis=0) * np.vstack([m.coef_[0] for m in cv_results["estimator"]])
        importances = StandardScaler(with_mean=False).fit(X).scale_ * model.coef_[0]
        df_importances = pd.DataFrame(importances, index=feature_names, columns=["importances"]).sort_values(
            "importances", ascending=False
        )

        # coef pvalues from formula
        # df_importances["pvalues"] = logit_pvalue(model, X.toarray())[1:]
    else:
        raise ValueError("Importance mode must be 'permutation' or 'coef'")

    return df_permutations, df_importances

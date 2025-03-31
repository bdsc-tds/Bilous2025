import sklearn
import numpy as np
import pandas as pd
import scipy
import gseapy
import anndata
from sklearn.model_selection import train_test_split, permutation_test_score
from sklearn.linear_model import LogisticRegression
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import StandardScaler
from typing import Dict, List

xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]


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
        if n_neighbors is not None:
            knn_graph = nn.kneighbors_graph()
        elif radius is not None:
            knn_graph = nn.radius_neighbors_graph()

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

    if np.sum(rnk > 0.0) == 0:
        print("Warning: all zero or all negative feature score vector. Returning empty dataframe.")
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


def get_mean_cell_identity_score_markers(
    ads: Dict[str, Dict[str, anndata.AnnData]],
    df_markers: pd.DataFrame,
    correction_methods: List[str],
    labels_key: str,
    columns=None,
) -> pd.DataFrame:
    """
    Calculates the mean number of genes (`n_genes`) for all unique cell types (`cti`) found
    in the `labels_key` column across different AnnData objects and correction methods.

    Args:
        ads (Dict[str, Dict[str, anndata.AnnData]]): A nested dictionary structure containing
            AnnData objects. The outer dictionary's keys are correction methods, and the
            inner dictionary's keys are sample identifiers. Each inner dictionary value
            is an AnnData object. Assumes `ads[correction_method][k]` is an AnnData object.
        df_markers (pd.DataFrame): A DataFrame containing cell type markers.
        correction_methods (List[str]): A list of correction methods to iterate through.
        labels_key (str): The key in `ad.obs` that contains cell type labels.
        columns (list): optional column names for the returned DataFrame
    Returns:
        pd.DataFrame: A DataFrame where each row represents a unique combination of sample
            identifier, correction method, and cell type. The columns include sample identifiers,
            correction method, cell type (`cti`), and the calculated mean number of genes
            ('n_genes').
    """

    data = []
    for correction_method in correction_methods:
        for k, ad in ads[correction_method].items():
            if ad is not None:
                unique_ctis = ad.obs[labels_key].unique()

                for cti in unique_ctis:
                    cti_markers = df_markers.query("cell_type == @cti")["gene_name"].tolist()
                    cti_markers = [g for g in cti_markers if g in ad.var_names]
                    if len(cti_markers):
                        n_genes_cti_markers = (ad[ad.obs[labels_key] == cti, cti_markers].X > 0).sum(1).A1
                        mean_n_genes = np.mean(n_genes_cti_markers)
                        data.append((*k, correction_method, cti, mean_n_genes))  # Append cell type to the data

    # Create the DataFrame from the collected data
    df = pd.DataFrame(data, columns=columns)

    if columns is not None:
        df.columns = columns
    return df


def get_mean_cell_identity_score_scrna(
    ads: Dict[str, Dict[str, anndata.AnnData]],
    ads_scrnaseq: Dict[str, anndata.AnnData],
    df_markers: pd.DataFrame,
    correction_methods: List[str],
    labels_key: str,
    columns=None,
) -> pd.DataFrame:
    """
    Calculates the mean number of genes (`n_genes`) for all unique cell types (`cti`) found
    in the `labels_key` column across different AnnData objects and correction methods.

    Args:
        ads (Dict[str, Dict[str, anndata.AnnData]]): A nested dictionary structure containing
            AnnData objects. The outer dictionary's keys are correction methods, and the
            inner dictionary's keys are sample identifiers. Each inner dictionary value
            is an AnnData object. Assumes `ads[correction_method][k]` is an AnnData object.
        ads_scrnaseq Dict[str, anndata.AnnData]: A dictionary structure containing
            AnnData objects. The dictionary's keys are tissue identifiers corresponding to keys in `ads`.
        correction_methods (List[str]): A list of correction methods to iterate through.
        labels_key (str): The key in `ad.obs` that contains cell type labels.
        columns (list): optional column names for the returned DataFrame
    Returns:
        pd.DataFrame: A DataFrame where each row represents a unique combination of sample
            identifier, correction method, and cell type. The columns include sample identifiers,
            correction method, cell type (`cti`), and the calculated mean number of genes
            ('n_genes').
    """

    data = []
    for correction_method in correction_methods:
        for k, ad in ads[correction_method].items():
            if ad is not None:
                unique_ctis = ad.obs[labels_key].unique()

                for cti in unique_ctis:
                    cti_markers = df_markers.query("cell_type == @cti")["gene_name"].tolist()
                    cti_markers = [g for g in cti_markers if g in ad.var_names]
                    if len(cti_markers):
                        n_genes_cti_markers = (ad[ad.obs[labels_key] == cti, cti_markers].X > 0).sum(1).A1
                        mean_n_genes = np.mean(n_genes_cti_markers)
                        data.append((*k, correction_method, cti, mean_n_genes))  # Append cell type to the data

    # Create the DataFrame from the collected data
    df = pd.DataFrame(data, columns=columns)

    if columns is not None:
        df.columns = columns
    return df


def get_cosine_similarity_score(pbs_xenium, pbs_scrna, labels_key, correction_methods, columns=None):
    """
    Calculates the cosine similarity score between each cell type's expression profile and its
    corresponding pseudobulk expression profile.

    Args:
        pbs_xenium (Dict[str, Dict[str, anndata.AnnData]]): A nested dictionary structure containing
            AnnData objects. The outer dictionary's keys are correction methods, and the
            inner dictionary's keys are sample identifiers. Each inner dictionary value
            is an AnnData object. Assumes `ads[correction_method][k]` is an AnnData object.
        pbs_scrna (Dict[str, anndata.AnnData]): A dictionary structure containing pseudobulk expression
            profiles. The dictionary's keys are tissue identifiers corresponding to keys in `pbs_xenium`.
        labels_key (str): The key in `pb_xenium.obs` that contains cell type labels.
        correction_methods (List[str]): A list of correction methods to iterate through.
        columns (list): optional column names for the returned DataFrame

    Returns:
        pd.DataFrame: A DataFrame where each row represents a unique combination of sample
            identifier, correction method, and cell type. The columns include sample identifiers,
            correction method, cell type (`cti`), and the calculated cosine similarity score
            ('cosine_sim').

    Notes:
        If a cell type is not found in the pseudobulk expression profile, a warning is printed.
    """
    data = []
    for correction_method in correction_methods:
        for k, pb_xenium in pbs_xenium[correction_method].items():
            print(correction_method, k)
            if pb_xenium is not None:
                tissue = k[1]
                pb_scrna = pbs_scrna[tissue]

                genes_tissue = np.intersect1d(pb_xenium.var_names, pb_scrna.var_names)
                unique_ctis = pb_xenium.obs_names.unique()
                for cti in unique_ctis:
                    if cti in pb_scrna.obs_names:
                        x = pb_xenium[pb_xenium.obs_names == cti, genes_tissue].X.toarray().squeeze()
                        y = pb_scrna[pb_scrna.obs_names == cti, genes_tissue].X.toarray().squeeze()

                        cosine_sim = 1 - scipy.spatial.distance.cosine(x, y)
                        data.append((*k, correction_method, cti, cosine_sim))  # Append cell type to the data
                    else:
                        print("Warning: cell type", cti, "not found in pseudobulk", tissue)

    # Create the DataFrame from the collected data
    df = pd.DataFrame(data, columns=columns)

    if columns is not None:
        df.columns = columns

    return df


def rename_correction_methods(df, column="correction_method"):
    df[column] = df[column].replace(
        {
            "resolvi": "ResolVI",
            "resolvi_supervised": "ResolVI supervised",
            "ovrlpy_correction_signal_integrity_threshold=0.5": "ovrlpy 0.5",
            "ovrlpy_correction_signal_integrity_threshold=0.7": "ovrlpy 0.7",
            "ovrlpy_0.5": "ovrlpy 0.5",
            "ovrlpy_0.7": "ovrlpy 0.7",
            "split_fully_purified": "SPLIT",
            "split": "SPLIT",
        }
    )


def extract_info_cell_type_pair(series, cti, ctj, flag):
    """Add info columns for number of cti cells being neighbor or not to ctj cells"""
    s = series.map(
        lambda df_: df_.query(f"label_key == '{cti}' and variable == 'has_{ctj}_neighbor' and value == {flag}")[
            "count"
        ].squeeze()
    )
    s = s.map(lambda x: 0 if isinstance(x, pd.Series) else x)  # replace empty series with 0

    return s


def get_df_summary_stats_plot(dfs, plot_metric="n_cells"):
    """
    Generate a summary DataFrame from a dictionary of summary statistics.

    Parameters:
    - dfs (dict): A dictionary containing summary statistics with correction methods as keys, as read by readwrite.read_diffexpr_results_samples.
    - plot_metric (str, optional): The metric to plot. Defaults to "n_cells". Can also be "df_has_neighbor_counts"
      or any other metric present in the summary statistics.

    Returns:
    - pd.DataFrame: A DataFrame containing the summary statistics, with columns based on the `xenium_levels`,
      correction method, and the specified plot metric. If the plot metric is a dictionary, additional columns are
      included for the dictionary keys and values.
    """

    data = []  # List to store data for the DataFrame

    for correction_method, k_values in dfs["summary_stats"].items():  # Iterate through correction methods
        for k, values in k_values.items():  # Iterate through keys k
            v = values[plot_metric]
            if plot_metric == "df_has_neighbor_counts":
                v = pd.DataFrame(v)
                row = list(k) + [correction_method, v]
                data.append(row)

            elif isinstance(v, dict):
                for k1, values1 in v.items():
                    row = list(k) + [correction_method, k1, values1]
                    data.append(row)

            else:
                row = list(k) + [correction_method, v]  # Create a row of data
                data.append(row)

    if isinstance(v, dict):
        columns = xenium_levels + ["correction_method", plot_metric + "_key", plot_metric + "_value"]
    else:
        columns = xenium_levels + ["correction_method", plot_metric]

    df = pd.DataFrame(data, columns=columns)
    rename_correction_methods(df)
    return df


def get_df_permutations_logreg_plot(df_permutations_logreg, correction_methods):
    df = {}
    for correction_method in correction_methods:
        for k, v in df_permutations_logreg[correction_method].items():
            for cti, ctj, _ in v.index:
                df[(correction_method, *k, cti, ctj)] = v.loc[(cti, ctj, 0)]
    df = pd.DataFrame(df).T.reset_index()
    df.columns = ["correction_method"] + xenium_levels + ["cti", "ctj"] + df.columns[-5:].tolist()
    rename_correction_methods(df)
    return df


def get_df_marker_rank_significance_plot(
    dfs_marker_rank_significance,
    rank_metric,
    plot_metric,
    correction_methods,
    use_precomputed,
):
    """
    Generate a DataFrame from a dictionary of DataFrames, to be used for plotting.

    Parameters:
    - dfs_marker_rank_significance (dict): A dictionary containing DataFrames with correction methods as keys, as read by readwrite.read_diffexpr_results_samples.
    - key (str): The key to use for accessing the DataFrames in dfs.
    - rank_metric (str): The ranking metric to use for the plot.
    - plot_metric (str): The metric to plot.
    - correction_methods (list): A list of correction methods to include in the plot.
    - use_precomputed (bool): Whether to use precomputed results (True) or not (False).

    Returns:
    - pd.DataFrame: A DataFrame containing the data for the plot, with columns for the correction method, xenium levels, cell types, and plot metric.
    """
    df = {}
    for correction_method in correction_methods:
        for k, v in dfs_marker_rank_significance[correction_method].items():
            if use_precomputed:
                rank_metric_ = rank_metric if correction_method == "raw" else rank_metric + "_precomputed"
            else:
                rank_metric_ = rank_metric
            df[(correction_method, *k)] = v.loc[rank_metric_, v.columns.get_level_values(2) == plot_metric]
    df = pd.concat(df).reset_index()
    df.columns = ["correction_method"] + xenium_levels + ["cti", "ctj", "plot_metric", plot_metric]
    rename_correction_methods(df)
    return df


# def pseudobulk(adata, labels_key):
#     """
#     Generate pseudo-bulk RNA-seq data from single-cell RNA-seq datasets.

#     Parameters
#     ----------
#     adata : Anndata
#     labels_key : str
#         The key in `adata.obs` used to identify cell type labels.

#     Returns
#     -------
#     pd.DataFrame
#         Contains pseudo-bulk expression profiles for each cell type,
#         with cell types as columns and genes as rows.
#     """

#     pbdata = {}

#     unique_cell_types = adata.obs[labels_key].unique()
#     for ct in unique_cell_types:
#         pb_cti_scrna = adata[adata.obs[labels_key] == ct].X.mean(0).A1
#         pbdata[ct] = pd.Series(pb_cti_scrna, index=adata.var_names)
#     pbdata = pd.DataFrame(pbdata)
#     return pbdata


def pseudobulk(ad, key, mode="sum"):
    """
    Generate pseudo-bulk RNA-seq data from single-cell RNA-seq datasets.

    Parameters
    ----------
    ad : Anndata
    key : str
        The key in `ad.obs` used to identify cell type labels.
    mode : str, optional
        Whether to sum or mean the expression values across cells of the same cell type.
        Options are "sum" (default) or "mean".

    Returns
    -------
    Anndata
        Contains pseudo-bulk expression profiles for each cell type,
        with cell types as columns and genes as rows.
    """

    def _agg_obs(x):
        val = x.mode()
        if len(val):
            return val[0]
        else:
            return None

    agg = {}
    for c in ad.obs[key].unique():
        if pd.isna(c):
            idx = ad.obs[key].isna()
        else:
            idx = ad.obs[key] == c
        if mode == "sum":
            agg[c] = np.asarray(ad[idx].X.sum(0)).squeeze()
        elif mode == "mean":
            agg[c] = np.asarray(ad[idx].X.mean(0)).squeeze()
    ad_states = anndata.AnnData(pd.DataFrame(agg).T)
    ad_states.var_names = ad.var_names
    ad_states.obs = ad.obs.groupby(key).agg(_agg_obs)
    return ad_states

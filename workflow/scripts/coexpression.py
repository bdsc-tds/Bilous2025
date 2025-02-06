import numpy as np
import pandas as pd
import scipy


def get_exclusive_gene_pairs(df_markers):
    """
    Get exclusive gene pairs from a DataFrame of cell type markers.
    Markers for a given cell type are assumed to be exclusive with all other cell types

    Parameters
    ----------
    df_markers : pd.DataFrame
        A DataFrame with first column corresponding to cell types and second column for gene names.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns 'cti', 'ctj', 'genei', 'genej'.
    """
    u_cell_types = df_markers.iloc[:, 0].unique()
    n_cell_types = len(u_cell_types)

    exclusive_gene_pairs = []
    for i in range(n_cell_types):
        for j in range(i + 1, n_cell_types):
            cti = u_cell_types[i]
            ctj = u_cell_types[j]

            df_markers_cti = df_markers[df_markers.iloc[:, 0] == cti]
            df_markers_ctj = df_markers[df_markers.iloc[:, 0] == ctj]

            for genei in df_markers_cti.iloc[:, 1]:
                for genej in df_markers_ctj.iloc[:, 1]:
                    exclusive_gene_pairs.append(
                        (*sorted((cti, ctj)), *sorted((genei, genej)))
                    )

    exclusive_gene_pairs = pd.DataFrame(
        exclusive_gene_pairs, columns=["cti", "ctj", "genei", "genej"]
    )
    return exclusive_gene_pairs


def counts_to_positivity(X, cutoff):
    """
    Convert counts to a binary positivity matrix based on a cutoff value.

    Parameters:
    X (scipy.sparse matrix): Input matrix with counts.
    cutoff (float): Cutoff value to determine positivity.

    Returns:
    tuple: A tuple containing the binary positivity matrix and positivity rate.
    """
    positivity = (X >= cutoff).astype(float)
    positivity_rate = positivity.mean(axis=0).A1
    # positivity_rate = np.clip(positivity.mean(axis=0).A1, min_pos_rate, None)
    return positivity, positivity_rate


def conditional_coexpression(positivity):
    """
    Compute conditional co-expression from a positivity matrix.

    Parameters:
    positivity (numpy.ndarray): Binary positivity matrix.

    Returns:
    numpy.ndarray: Conditional co-expression matrix.
    """
    coex = positivity.T @ positivity
    cond_coex = coex / positivity.sum(axis=0)
    return cond_coex.toarray()


def jaccard_coexpression(X):
    """
    Compute the Jaccard index between the rows of X.
    Parameters:
    X (scipy.sparse matrix): Input binary matrix.

    Returns:
    numpy.ndarray: Jaccard index matrix.
    """

    X = X.astype(bool).astype(int)
    intrsct = X.dot(X.T)
    row_sums = intrsct.diagonal()
    unions = row_sums[:, None] + row_sums - intrsct
    jaccard_index = intrsct / unions

    # if min_jaccard > 0.0:
    #     min_jaccard = np.nan_to_num(min_jaccard)
    #     min_jaccard = np.clip(jaccard_index, min_jaccard, None)
    return jaccard_index


def pearson_coexpression(X):
    """
    Compute the Pearson correlation coefficient matrix.

    Parameters:
    X (scipy.sparse matrix): Input matrix.

    Returns:
    numpy.ndarray: Pearson correlation coefficient matrix.
    """
    return np.corrcoef(X.toarray().T)


def spearman_coexpression(X):
    """
    Compute the Spearman rank correlation coefficient matrix.

    Parameters:
    X (scipy.sparse matrix): Input matrix.

    Returns:
    numpy.ndarray: Spearman correlation coefficient matrix.
    """
    return scipy.stats.spearmanr(X.toarray()).statistic


# def odds_ratio_coexpression(positivity, min_odds_ratio=0.01):
#     # Compute joint and marginal sums
#     joint_counts = (
#         positivity.T @ positivity
#     )  # Joint presence counts for each pair of genes
#     total_cells = positivity.shape[0]

#     gene_sums = positivity.sum(axis=0)
#     neither_counts = (
#         total_cells - gene_sums[:, None] - gene_sums + joint_counts
#     )  # Cells with neither gene expressed
#     gene_only_counts_a = (
#         gene_sums[:, None] - joint_counts
#     )  # Cells with only gene A expressed
#     gene_only_counts_b = gene_sums - joint_counts  # Cells with only gene B expressed

#     # Calculate odds ratio matrix
#     odds_ratio_matrix = (joint_counts * neither_counts) / (
#         gene_only_counts_a * gene_only_counts_b + 1e-10
#     )

#     # Clip odds ratios to avoid values below the minimum threshold
#     return np.clip(odds_ratio_matrix, min_odds_ratio, np.inf)


# def relative_risk_coexpression(positivity, min_relative_risk=0.01):
#     # Compute joint and marginal sums
#     joint_counts = (
#         positivity.T @ positivity
#     )  # Joint presence counts for each pair of genes
#     total_cells = positivity.shape[0]

#     gene_sums = positivity.sum(axis=0)

#     # Calculate probabilities
#     prob_joint = joint_counts / total_cells  # Probability of both genes expressed
#     prob_a = gene_sums / total_cells  # Probability of gene A expressed
#     prob_b = gene_sums / total_cells  # Probability of gene B expressed

#     # Calculate relative risk matrix
#     relative_risk_matrix = prob_joint / (prob_a * prob_b + 1e-10)

#     # Clip relative risks to avoid values below the minimum threshold
#     return np.clip(relative_risk_matrix, min_relative_risk, np.inf)


# def relative_coexpression(CC_ref_seg, cc, gene_pairs):
#     rel_ccs = []
#     for a, b in gene_pairs:
#         i = CC_ref_seg.genes.index(a)
#         j = CC_ref_seg.genes.index(b)
#         CC_ref_seg_ij = CC_ref_seg.CC[i, j]

#         i = cc.genes.index(a)
#         j = cc.genes.index(b)
#         cc_ij = cc.CC[i, j]

#         rel_cc = np.log2(cc_ij / CC_ref_seg_ij)
#         rel_ccs.append(rel_cc)

#     return rel_ccs


# def thin_counts_v2(X, target_size, rs=None):
#     if rs is None:
#         rs = np.random.RandomState()

#     n_counts = X.sum(axis=1)
#     probabilities = (X / n_counts).tocsr()

#     X_downsampled = X.copy()
#     for i in range(X.shape[0]):
#         diff_size = n_counts[i] - target_size
#         if diff_size > 0:
#             X_downsampled[i] -= rs.multinomial(
#                 n=diff_size, pvals=probabilities[i].toarray().flat
#             )

#     # Ensure no negative values by recursion if needed
#     if np.any(X_downsampled.data < 0):
#         print("clipping negative counts to 0...")
#         X_downsampled.data = np.clip(X_downsampled.data, 0, None, out=X_downsampled)
#         return thin_counts(X_downsampled, target_size=target_size)

#     return X_downsampled


# def thin_counts_v3(X, target_size, rs=None):
#     if rs is None:
#         rs = np.random.RandomState()

#     n_counts = X.sum(axis=1)
#     probabilities = (X / n_counts).tocsr()

#     X_downsampled = X.copy()
#     for i in tqdm(range(X.shape[0])):
#         diff_size = n_counts[i] - target_size
#         pvals = probabilities[i].toarray().flat
#         while diff_size > 0:
#             thin = rs.multinomial(n=1, pvals=pvals)
#             if X_downsampled[i, np.where(thin)[0]] > 0:
#                 X_downsampled[i] -= thin
#                 diff_size -= 1

#     return X_downsampled


# def thin_counts_old(X, target_count, rs=None):
#     if rs is None:
#         rs = np.random.RandomState()

#     n_counts = X.sum(axis=1)
#     probabilities = (X / n_counts).tocsr()
#     rows = []
#     for i in range(probabilities.shape[0]):
#         counts = rs.multinomial(target_count, probabilities[i].toarray().flat)
#         rows.append(scipy.sparse.csr_matrix(counts))
#     return scipy.sparse.vstack(rows)


def thin_counts(X, target_count, gen=None):
    """
    Downdonor counts to a target count per row.

    Parameters:
    X (scipy.sparse matrix): Input count matrix.
    target_count (int): Target count for each row.
    gen (numpy.random.Generator, optional): Random generator instance.

    Returns:
    scipy.sparse.csr_matrix: Downdonord count matrix.
    """
    if gen is None:
        gen = np.random.default_rng()

    n_counts = X.sum(axis=1)
    probabilities = (X / n_counts).toarray()
    X_thin = np.random.default_rng().multinomial(target_count, probabilities)
    return scipy.sparse.csr_matrix(X_thin)


def sparsify(X):
    """
    Convert a dense matrix to a sparse matrix if not already sparse.

    Parameters:
    X (numpy.ndarray or scipy.sparse matrix): Input matrix.

    Returns:
    scipy.sparse.csr_matrix: Sparse matrix.
    """
    if not scipy.sparse.issparse(X):
        return scipy.sparse.csr_matrix(X.astype(float))
    return X.astype(float)


def coexpression(
    adata,
    positivity_cutoff=1,
    min_donors=0,
    target_count=50,
    method="conditional",
    seed=0,
    # min_positivity_rate: float = 0.01,
    # min_cond_coex: float = 0.0,
):
    """
    Calculate co-expression matrix using different methods.

    Parameters:
    adata (anndata.AnnData): AnnData object containing gene expression data.
    positivity_cutoff (float): Cutoff for determining positivity.
    min_donors (int): Minimum number of donors required.
    target_count (int): Target count for downsampling.
    method (str): Method for co-expression calculation.
    seed (int): Seed for random number generation.

    Returns:
    tuple: A tuple containing co-expression matrix, downdonord matrix, positivity matrix, positivity rate, and mask.
    """
    gen = np.random.default_rng(seed)

    X = sparsify(adata.X)

    # Apply mask based on target count threshold
    n_counts = X.sum(axis=1).A1

    if target_count is not None:
        mask = n_counts >= target_count
        print(
            sum(mask),
            "/",
            X.shape[0],
            "(",
            round(100 * sum(mask) / X.shape[0], 2),
            "% ) cells reaching the target count",
        )

        # Modify mask based on min_donors threshold
        if sum(mask) < min_donors:
            print(
                f"Less than {min_donors=} reach the target count. Setting to {min_donors}"
            )
            mask = np.argsort(n_counts)[::-1][:min_donors]

        # Downdonor counts
        X_downsample = thin_counts(X[mask], target_count, gen=gen)
    else:
        mask = np.ones(X.shape[0], dtype=bool)
        X_downsample = X

    # Convert counts to binary positivity matrix based on threshold
    pos, pos_rate = counts_to_positivity(
        X_downsample,
        positivity_cutoff,
    )  # min_positivity_rate)

    # Calculate conditional co-expression
    if method == "conditional":
        CC = conditional_coexpression(pos)
    elif method == "jaccard":
        CC = jaccard_coexpression(pos.T).toarray()
    elif method == "pearson":
        CC = pearson_coexpression(X_downsample)
    elif method == "spearman":
        CC = spearman_coexpression(X_downsample)

    CC = pd.DataFrame(CC, index=adata.var_names, columns=adata.var_names)
    pos_rate = pd.Series(pos_rate, index=adata.var_names)

    return CC, X_downsample, pos, pos_rate, mask


def censored_ratio(
    CC_ref_seg,
    CC_other_seg,
    pos_rate_ref_seg=None,
    pos_rate_other_seg=None,
    min_positivity_rate=0.0,
    min_cond_coex=0.0,
    log2=True,
):
    """
    Compute the ratio of co-expression matrices with optional censoring and log transformation.

    Parameters:
    CC_ref_seg (pd.DataFrame): Reference co-expression matrix.
    CC_other_seg (pd.DataFrame): Other co-expression matrix for comparison.
    pos_rate_ref_seg (pd.Series, optional): Positivity rate for reference segmentation.
    pos_rate_other_seg (pd.Series, optional): Positivity rate for other segmentation.
    min_positivity_rate (float): Minimum positivity rate for filtering.
    log2 (bool): Whether to apply log2 transformation.

    Returns:
    pd.DataFrame: Censored and transformed ratio matrix.
    """
    # pseudocount to avoid zero numerators or denominators
    if min_cond_coex > 0.0:
        CCdiff = CC_other_seg.clip(lower=min_cond_coex) / CC_ref_seg.clip(
            lower=min_cond_coex
        )
    else:
        CCdiff = CC_other_seg / CC_ref_seg

    if log2:
        CCdiff = np.log2(CCdiff)

    # exclude stuff that's barely expressed in one or the other
    if min_positivity_rate > 0.0:
        mask = (pos_rate_ref_seg < min_positivity_rate) | (
            pos_rate_other_seg < min_positivity_rate
        )
        CCdiff.loc[mask, mask] = 0.0

    return CCdiff


def compare_segmentations(
    CC_ref_seg,
    CC_other_seg,
    pos_rate_ref_seg,
    pos_rate_other_seg,
    min_positivity_rate=0.01,
    min_cond_coex=0.05,
    cc_cutoff=1.5,
    method=None,
    log2=True,
):
    """
    Compare two co-expression segmentations and identify spurious gene pairs.
    Parameters:
    CC_ref_seg (pd.DataFrame): Reference co-expression matrix.
    CC_other_seg (pd.DataFrame): Other co-expression matrix for comparison.
    pos_rate_ref_seg (pd.Series): Positivity rate for reference segmentation.
    pos_rate_other_seg (pd.Series): Positivity rate for other segmentation.
    min_positivity_rate (float): Minimum positivity rate for filtering.
    cc_cutoff (float): Cutoff for spurious gene pair identification.
    method (str, optional): Method for co-expression calculation.
    log2 (bool): Whether to apply log2 transformation.

    Returns:
    tuple: A tuple containing the difference matrix and spurious gene pairs.
    """

    CCdiff = censored_ratio(
        CC_ref_seg,
        CC_other_seg,
        pos_rate_ref_seg=pos_rate_ref_seg,
        pos_rate_other_seg=pos_rate_other_seg,
        min_positivity_rate=min_positivity_rate,
        min_cond_coex=min_cond_coex,
        log2=log2,
    )

    if log2:
        cc_cutoff = np.log2(cc_cutoff)

    if method == "conditional":
        spurious_gene_pairs = np.where(CCdiff >= cc_cutoff)
    else:
        CCdiff_triu = np.triu(CCdiff, 1)
        spurious_gene_pairs = np.where(CCdiff_triu >= cc_cutoff)

    spurious_gene_pairs = np.array(
        (
            CCdiff.index[spurious_gene_pairs[0]],
            CCdiff.columns[spurious_gene_pairs[1]],
        )
    ).T
    return CCdiff, spurious_gene_pairs


def coexpression_by_cell_type(CC, genes_markers, df_markers_panel):
    """
    Calculate co-expression by cell type based on marker genes.

    Parameters:
    CC (pd.DataFrame): Co-expression matrix.
    genes_markers (list): List of marker genes.
    df_markers_panel (pd.DataFrame): DataFrame of markers and corresponding cell types.

    Returns:
    tuple: DataFrames of co-expression by cell type and melted co-expression.
    """
    df_cc_melt = (
        CC.loc[genes_markers, genes_markers].melt(ignore_index=False).reset_index()
    )
    df_cc_melt = df_cc_melt[df_cc_melt["index"] != df_cc_melt["variable"].values]
    df_cc_melt.columns = ["genei", "genej", "value"]

    u_ct = df_markers_panel.iloc[:, 0].unique()
    df_markers_panel_ct = df_markers_panel.set_index("canonical").iloc[:, 0]
    df_ct_markers_cc = pd.DataFrame(0.0, index=u_ct, columns=u_ct)
    df_ct_markers_cc_count = pd.DataFrame(0.0, index=u_ct, columns=u_ct)

    for i, row in df_cc_melt.iterrows():
        ct_genei = df_markers_panel_ct[row["genei"]]
        ct_genej = df_markers_panel_ct[row["genej"]]
        df_ct_markers_cc.loc[ct_genei, ct_genej] += row["value"]
        df_ct_markers_cc_count.loc[ct_genei, ct_genej] += 1

    df_ct_markers_cc /= df_ct_markers_cc_count

    return df_ct_markers_cc, df_cc_melt


def coexpression_cells(pos_nuc):
    """
    Calculate co-expression for each cell from a positivity matrix.

    Parameters:
    pos_nuc (numpy.ndarray): Positivity matrix for nuclei.

    Returns:
    list: List of co-expression matrices for each cell.
    """
    CC_cells = []
    for i in range(pos_nuc.shape[0]):
        CC_cells.append(pos_nuc[i].T @ pos_nuc[i])
    return CC_cells


def coexpression_cells_score(CC_cells, marker_genes_idx):
    """
    Calculate co-expression scores for each cell based on marker genes.

    Parameters:
    CC_cells (list): List of co-expression matrices for each cell.
    marker_genes_idx (list): List of indices for marker genes.

    Returns:
    numpy.ndarray: Array of co-expression scores for each cell.
    """
    CC_cells_score = []
    for i in range(len(CC_cells)):
        CC_cells_score.append(CC_cells[i][marker_genes_idx][:, marker_genes_idx].sum())

    return np.array(CC_cells_score)


def coexpression_cells_score_gene_pairs(CC_cells, gene_pairs_idx):
    """
    Calculate co-expression scores for each cell based on gene pairs.

    Parameters:
    CC_cells (list): List of co-expression matrices for each cell.
    gene_pairs_idx (numpy.ndarray): Array of indices for gene pairs.

    Returns:
    numpy.ndarray: Array of co-expression scores for each cell.
    """
    CC_cells_score = []
    for i in range(len(CC_cells)):
        CC_cells_score.append(
            CC_cells[i][gene_pairs_idx[:, 0], gene_pairs_idx[:, 1]].sum()
        )
    return np.array(CC_cells_score)


def find_markers(adata, ct_key, threshold_fraction=0.05, threshold_diff=0.25):
    """
    Identify markers from an anndata based on mutually exclusive gene expression in each cell type

    Parameters:
        adata (anndata.AnnData):
            AnnData with cell type labels.
        ct_key (str):
            adata.obs column name for cell type labels.
        threshold_fraction (float):
            Maximum fraction of expressing cells allowed in negative cell types.
        threshold_diff (float):
            Minimum difference in fraction of expressing cells
            to consider between cell type of interest and negative cell type.

    Returns:
        df_max_difference_thresholded (pd.DataFrame):
            Max difference (best marker gene score) found for each cell type pair.
        df_difference_thresholded (pd.DataFrame):
            Scores for each cell type pair and each gene
    """
    u_cell_types = adata.obs[ct_key].unique()

    # get df with fraction of cells expressing each gene for each cell type
    df_fraction = pd.DataFrame(index=adata.var_names, columns=u_cell_types)
    for ct in u_cell_types:
        adata_ct = adata[adata.obs[ct_key] == ct]
        df_fraction[ct] = (adata_ct.X > 0).mean(0).A1

    # Initialize result DataFrame and dictionary
    df_max_difference_thresholded = pd.DataFrame(
        0.0, index=u_cell_types, columns=u_cell_types
    )
    df_difference_thresholded = {}

    # Compute differences and thresholds
    for cti in u_cell_types:
        for ctj in u_cell_types:
            if cti == ctj:
                continue
            else:
                diff = (df_fraction[cti] - df_fraction[ctj]).clip(0)
                diff[df_fraction[ctj] > threshold_fraction] = 0.0
                diff[diff < threshold_diff] = 0.0

                df_max_difference_thresholded.loc[
                    ctj, cti
                ] = df_max_difference_thresholded.loc[cti, ctj] = diff.max()
                df_difference_thresholded[cti, ctj] = diff

    # Convert dictionary to DataFrame and sort columns by maximum values
    df_difference_thresholded = pd.DataFrame(df_difference_thresholded)
    df_difference_thresholded.columns = [
        "_".join(pair) for pair in df_difference_thresholded.columns
    ]

    return df_fraction, df_difference_thresholded, df_max_difference_thresholded

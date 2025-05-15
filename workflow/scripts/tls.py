import numpy as np
import pandas as pd
import anndata
from sklearn.neighbors import radius_neighbors_graph
from sklearn.cluster import HDBSCAN  # Import sklearn's HDBSCAN directly
from collections import Counter
import time


# --- Function 1: TLS Detection (Fixed % Thresholds, sklearn HDBSCAN) ---
def detect_tls(
    adata: anndata.AnnData,
    radius: float = 30.0,
    b_cell_label: str = "B cell",
    t_cell_label: str = "T cell",
    min_perc_b: float = 0.0,  # Added: Min % B cells in neighborhood (0.0-1.0)
    min_perc_t: float = 0.0,  # Added: Min % T cells in neighborhood (0.0-1.0)
    fold_change_threshold: float = 1.5,
    min_b: int = 5,  # Min absolute number of B cells in neighborhood
    min_t: int = 5,  # Min absolute number of T cells in neighborhood
    min_tls_size: int = 20,  # Used for HDBSCAN min_cluster_size AND post-filter B+T count
    cell_type_obs_key: str = "cell_type",
    spatial_obsm_key: str = "spatial",
    output_obs_key: str = "tls_cluster",
    n_jobs: int = 1,
) -> anndata.AnnData:
    """
    Detects potential TLS using fixed percentage and absolute count thresholds
    for neighborhood selection, followed by sklearn's HDBSCAN clustering and
    a post-clustering B+T cell count filter.

    Method:
    1. Build radius neighbors graph.
    2. Calculate B/T cell counts and proportions in each neighborhood.
    3. Identify 'candidate' cells whose neighborhoods meet minimum percentage AND
       minimum absolute count thresholds for both B and T cells.
    4. Cluster spatial coordinates of candidate cells using sklearn.cluster.HDBSCAN,
       using `min_tls_size` for the `min_cluster_size` parameter.
    5. Filter resulting HDBSCAN clusters based on `min_tls_size` (total B+T cells).
    6. Annotate valid TLS clusters in adata.obs.

    Args:
        adata: Annotated data matrix.
        radius: Radius for neighborhood definition.
        b_cell_label, t_cell_label: Labels for B/T cells.
        min_perc_b: Minimum proportion (0.0-1.0) of B cells required in neighborhood.
        min_perc_t: Minimum proportion (0.0-1.0) of T cells required in neighborhood.
        min_b: Minimum absolute number of B cells required in neighborhood.
        min_t: Minimum absolute number of T cells required in neighborhood.
        min_tls_size: Used for HDBSCAN `min_cluster_size` AND the minimum total B+T cells
                      required within a final cluster post-filtering.
        cell_type_obs_key, spatial_obsm_key, output_obs_key: AnnData keys.
        n_jobs: Number of parallel jobs for radius_neighbors_graph.
        hdbscan_min_samples: Optional `min_samples` parameter for HDBSCAN. If None,
                             defaults based on min_cluster_size internally.

    Returns:
        AnnData object with `adata.obs[output_obs_key]`.
    """
    print("--- Starting TLS Detection (Fixed % / Abs Count Thresholds, sklearn HDBSCAN) ---")
    start_time_local = time.time()

    # --- Input Validation ---
    if spatial_obsm_key not in adata.obsm:
        raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_obsm_key}']")
    if cell_type_obs_key not in adata.obs:
        raise ValueError(f"Cell type annotations not found in adata.obs['{cell_type_obs_key}']")

    coords = adata.obsm[spatial_obsm_key].astype(np.float64)
    cell_types = adata.obs[cell_type_obs_key]

    # Update parameter print statement
    print(
        f"Params: {radius=}, "
        f"min_perc_b={min_perc_b:.2f}, min_perc_t={min_perc_t:.2f}, "
        f"min_abs_count_b={min_b}, min_abs_count_t={min_t}, "
        f"{min_tls_size=} (used for HDBSCAN min_cluster_size & final filter)"
    )

    # --- Global Proportions (Calculated but not used for filtering) ---
    global_proportions = cell_types.value_counts(normalize=True)
    global_prop_b = global_proportions.get(b_cell_label, 0.0)
    global_prop_t = global_proportions.get(t_cell_label, 0.0)
    print(f"Global props (for info): B={global_prop_b:.4f}, T={global_prop_t:.4f}")

    # --- Radius Neighbors Graph ---
    print(f"Building detection radius neighbors graph (radius={radius}, n_jobs={n_jobs})...")
    adjacency_graph = radius_neighbors_graph(
        coords, radius=radius, mode="connectivity", include_self=True, n_jobs=n_jobs
    )

    # --- Neighborhood Composition ---
    print("Calculating detection neighborhood compositions...")
    is_b_cell = (cell_types == b_cell_label).values
    is_t_cell = (cell_types == t_cell_label).values
    neigh_count_b = adjacency_graph @ is_b_cell.astype(int)
    neigh_count_t = adjacency_graph @ is_t_cell.astype(int)
    total_neighbors = np.array(adjacency_graph.sum(axis=1)).flatten()

    # Calculate proportions safely
    neigh_prop_b = np.zeros_like(neigh_count_b, dtype=float)
    neigh_prop_t = np.zeros_like(neigh_count_t, dtype=float)
    valid_neighbors_mask = total_neighbors > 0
    neigh_prop_b[valid_neighbors_mask] = np.divide(
        neigh_count_b[valid_neighbors_mask], total_neighbors[valid_neighbors_mask]
    )
    neigh_prop_t[valid_neighbors_mask] = np.divide(
        neigh_count_t[valid_neighbors_mask], total_neighbors[valid_neighbors_mask]
    )
    # Note: Cells with 0 neighbors will have prop 0.0, which is usually fine.

    # --- Identify Candidate Cells (Using Fixed % and Abs Count Thresholds) ---
    print("Identifying candidate cells based on fixed thresholds...")
    # Fixed percentage checks
    perc_b_check = neigh_prop_b >= min_perc_b
    perc_t_check = neigh_prop_t >= min_perc_t

    # Absolute *COUNT* checks
    absolute_b_check = neigh_count_b >= min_b
    absolute_t_check = neigh_count_t >= min_t

    # Combine all conditions
    adaptive_b_check = (global_prop_b == 0) | (neigh_prop_b > (fold_change_threshold * global_prop_b))
    adaptive_t_check = (global_prop_t == 0) | (neigh_prop_t > (fold_change_threshold * global_prop_t))
    absolute_b_check = neigh_count_b >= min_b
    absolute_t_check = neigh_count_t >= min_t
    is_candidate = (
        perc_b_check & absolute_b_check & adaptive_b_check & perc_t_check & absolute_t_check & adaptive_t_check
    )

    adata.obs["neigh_count_b"] = neigh_count_b  # Keep for potential inspection
    adata.obs["neigh_count_t"] = neigh_count_t
    adata.obs["neigh_prop_b"] = neigh_prop_b
    adata.obs["neigh_prop_t"] = neigh_prop_t
    adata.obs["is_candidate"] = is_candidate

    candidate_indices = np.where(is_candidate)[0]
    n_candidates = len(candidate_indices)
    print(f"Found {n_candidates} candidate cells.")

    # --- Cluster Candidates (sklearn HDBSCAN) ---
    adata.obs[output_obs_key] = -1  # Initialize output
    # Compare n_candidates against the clustering size parameter
    if n_candidates < min_tls_size:
        print(f"Skipping clustering: candidates ({n_candidates}) < HDBSCAN min_cluster_size ({min_tls_size}).")
        elapsed = time.time() - start_time_local
        print(f"--- TLS Detection Finished ({elapsed:.2f} seconds) ---")
        return adata

    print(f"Clustering candidate cells using sklearn.cluster.HDBSCAN (min_cluster_size={min_tls_size})...")
    candidate_coords = coords[candidate_indices]

    try:
        # Use sklearn's HDBSCAN
        hdb = HDBSCAN(
            min_cluster_size=min_tls_size,
            metric="euclidean",
        )
        cluster_labels = hdb.fit_predict(candidate_coords)
    except Exception as e:  # Catch clustering errors
        print(f"Error during HDBSCAN clustering: {e}")
        print("Clustering failed. No TLS identified.")
        elapsed = time.time() - start_time_local
        print(f"--- TLS Detection Finished ({elapsed:.2f} seconds) ---")
        return adata

    n_clusters_found = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise_points = np.sum(cluster_labels == -1)
    print(f"HDBSCAN found {n_clusters_found} initial clusters and {n_noise_points} noise points.")

    if n_clusters_found == 0:
        print("No dense clusters identified by HDBSCAN.")
        elapsed = time.time() - start_time_local
        print(f"--- TLS Detection Finished ({elapsed:.2f} seconds) ---")
        return adata

    # --- Filter Clusters by Size (Post-Clustering B+T count check using min_tls_size) ---
    # This filter remains as per the original code block provided, using the same threshold
    print(f"Filtering clusters by minimum biological size (>= {min_tls_size} B+T cells)...")
    tls_counter = 0
    unique_cluster_ids = sorted([label for label in np.unique(cluster_labels) if label != -1])

    # Use original adata index for mapping
    original_indices = adata.obs.index[candidate_indices]
    cluster_map = pd.Series(cluster_labels, index=original_indices)

    for cluster_id in unique_cluster_ids:
        cluster_member_indices_original = cluster_map[cluster_map == cluster_id].index
        cluster_cell_types = cell_types.loc[cluster_member_indices_original]

        cluster_counts = Counter(cluster_cell_types)
        num_b_in_cluster = cluster_counts.get(b_cell_label, 0)
        num_t_in_cluster = cluster_counts.get(t_cell_label, 0)
        total_bt_in_cluster = num_b_in_cluster + num_t_in_cluster

        # Apply the biological size filter (`min_tls_size`)
        if total_bt_in_cluster >= min_tls_size:
            # Assign the new sequential label using original index labels
            adata.obs.loc[cluster_member_indices_original, output_obs_key] = tls_counter
            print(f"  Cluster {cluster_id} PASSED: {total_bt_in_cluster} B+T cells. Assigned TLS label {tls_counter}.")
            tls_counter += 1
        # else: # Optional verbose output for failed clusters
        #     print(f"  Cluster {cluster_id} FAILED: {total_bt_in_cluster} B+T cells (min_tls_size = {min_tls_size}).")

    print(f"Identified {tls_counter} potential TLS regions after size filtering.")
    elapsed = time.time() - start_time_local
    print(f"--- TLS Detection Finished ({elapsed:.2f} seconds) ---")

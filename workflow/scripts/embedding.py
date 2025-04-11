import importlib
import scvi
import scanpy as sc

if importlib.util.find_spec("rapids_singlecell") is not None:
    # Import gpu libraries, Initialize rmm and cupy
    import rapids_singlecell as rsc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator

    rmm.disable_logging()


def check_gpu_availability():
    import torch


def umap(
    adata,
    use_rep="X_pca",
    n_neighbors=15,
    metric="euclidean",
    min_dist=0.5,
    backend="gpu",
):
    neighbors_params_keys = [
        "n_neighbors",
        "method",
        "random_state",
        "metric",
        "use_rep",
    ]

    neighbors_params = dict(zip(neighbors_params_keys, [n_neighbors, "umap", 0, metric, use_rep]))

    if "neighbors" not in adata.uns:
        compute_neighbors = True
    else:
        compute_neighbors = any(
            [adata.uns["neighbors"]["params"][k] != neighbors_params[k] for k in neighbors_params_keys]
        )

    if backend == "gpu" and check_gpu_availability():
        print("Transferring data to GPU...")
        rsc.get.anndata_to_GPU(adata)

        if compute_neighbors:
            rsc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, metric=metric)
        adata.obsm[f"{use_rep}_umap"] = rsc.tl.umap(adata, copy=True, min_dist=min_dist).obsm["X_umap"]

        print("Transferring data back to CPU...")
        rsc.get.anndata_to_CPU(adata)
    else:
        if backend == "gpu":
            print("GPU not available. Switching to CPU backend...")

        if compute_neighbors:
            sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, metric=metric)
        sc.tl.umap(adata, min_dist=min_dist)
        adata.obsm[f"{use_rep}_umap"] = adata.obsm["X_umap"]


def mde(
    adata,
    use_rep="X_pca",
):
    adata.obsm[f"{use_rep}_mde"] = scvi.model.utils.mde(adata.obsm[use_rep])

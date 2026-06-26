import anndata as ad
import scipy.io
import scipy.sparse as sp
import pandas as pd
import numpy as np
import os

# ── USER CONFIG ───────────────────────────────────────────────────────────────
H5AD_PATH = "/path/to/external_lung.h5ad"
OUT_DIR   = "/path/to/output/external_lung_mtx"
# ─────────────────────────────────────────────────────────────────────────────

adata = ad.read_h5ad(H5AD_PATH)

os.makedirs(OUT_DIR, exist_ok=True)

# Raw counts stored in 'count' layer; fall back to X
mat = adata.layers["count"] if "count" in adata.layers else adata.X
if not sp.issparse(mat):
    mat = sp.csr_matrix(mat)

scipy.io.mmwrite(f"{OUT_DIR}/matrix.mtx", mat.T)

# Gene symbols live in var['feature_name'] (stored as categorical)
gene_symbols = adata.var["feature_name"].values
pd.Series(adata.obs_names).to_csv(f"{OUT_DIR}/barcodes.txt", index=False, header=False)
pd.Series(gene_symbols).to_csv(f"{OUT_DIR}/features.txt", index=False, header=False)

adata.obs.to_csv(f"{OUT_DIR}/metadata.csv")

if "X_umap" in adata.obsm:
    pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"]
    ).to_csv(f"{OUT_DIR}/umap.csv")

print("Done — files written to", OUT_DIR)

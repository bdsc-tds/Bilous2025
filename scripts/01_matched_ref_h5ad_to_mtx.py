import anndata as ad
import scipy.io
import scipy.sparse as sp
import pandas as pd
import os

# ── USER CONFIG ───────────────────────────────────────────────────────────────
H5AD_PATH = "/path/to/lung.h5ad"
OUT_DIR   = "/path/to/output/lung_mtx"
# ─────────────────────────────────────────────────────────────────────────────

adata = ad.read_h5ad(H5AD_PATH)

out = OUT_DIR
os.makedirs(out, exist_ok=True)

mat = adata.layers["SoupX"] if "SoupX" in adata.layers else adata.X
if not sp.issparse(mat):
    mat = sp.csr_matrix(mat)

scipy.io.mmwrite(f"{out}/matrix.mtx", mat.T)

# Plain barcode / feature lists (no header) for ReadMtx
pd.Series(adata.obs_names).to_csv(f"{out}/barcodes.txt", index=False, header=False)
pd.Series(adata.var_names).to_csv(f"{out}/features.txt", index=False, header=False)

# Full metadata separately
adata.obs.to_csv(f"{out}/metadata.csv")

if "X_umap" in adata.obsm:
    pd.DataFrame(
        adata.obsm["X_umap"],
        index=adata.obs_names,
        columns=["UMAP_1", "UMAP_2"]
    ).to_csv(f"{out}/umap.csv")

print("Done — files written to", out)

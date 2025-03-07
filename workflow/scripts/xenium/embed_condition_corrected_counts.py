import dask

dask.config.set({"dataframe.query-planning": False})

from pathlib import Path
import argparse
import scipy
import pandas as pd
import scanpy as sc
import sys

sys.path.append("workflow/scripts/")
import preprocessing

# Set up argument parser
parser = argparse.ArgumentParser(description="Embed panel of Xenium donors.")
parser.add_argument("--condition", type=Path, help="Path to the panel file.")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--normalisation", type=str, help="Normalisation method")
parser.add_argument("--layer", type=str, help="Name of saved layer of the seurat object, data or scale_data")
parser.add_argument("--n_comps", type=int, help="Number of components.")
parser.add_argument("--n_neighbors", type=int, help="Number of neighbors.")
parser.add_argument("--metric", type=str, help="Distance metric to use.")
parser.add_argument("--min_dist", type=float, help="Minimum distance parameter.")
parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
parser.add_argument(
    "--raw_corrected_counts", action="store_true", help="load raw corrected counts .h5 instead of normalised counts"
)

args = parser.parse_args()

# Access the arguments
condition = args.condition
out_file = args.out_file
normalisation = args.normalisation
layer = args.layer
n_comps = args.n_comps
n_neighbors = args.n_neighbors
metric = args.metric
min_dist = args.min_dist
min_counts = args.min_counts
min_features = args.min_features
max_counts = args.max_counts
max_features = args.max_features
min_cells = args.min_cells
raw_corrected_counts = args.raw_corrected_counts

correction_method = condition.parents[1].stem
segmentation = condition.parents[0].stem

# read xenium samples
if raw_corrected_counts:
    ads = {}
    for panel in (panels := condition.iterdir()):
        for donor in (donors := panel.iterdir()):
            for sample in (samples := donor.iterdir()):
                k = (segmentation, condition, panel.stem, donor.stem, sample.stem)
                sample_path = sample / "corrected_counts.h5"
                ads[k] = sc.read_10x_h5(sample_path)
    is_raw = True

else:
    ads = {}
    for panel in (panels := condition.iterdir()):
        for donor in (donors := panel.iterdir()):
            for sample in (samples := donor.iterdir()):
                print(sample)

                k = (segmentation, condition.stem, panel.stem, donor.stem, sample.stem)
                sample_counts_path = sample / f"{normalisation}/normalised_counts/{layer}.parquet"
                sample_idx_path = sample / f"{normalisation}/normalised_counts/cells.parquet"

                ads[k] = sc.AnnData(pd.read_parquet(sample_counts_path))
                if layer != "scale_data":  # no need to sparsify scale_data which is dense
                    ads[k].X = scipy.sparse.csr_matrix(ads[k].X)
                ads[k].obs_names = pd.read_parquet(sample_idx_path).iloc[:, 0]
    is_raw = False
    min_counts = None
    min_genes = None
    max_counts = None
    max_genes = None
    min_cells = None

# concatenate
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
for k in ads.keys():
    for i, l in enumerate(xenium_levels):
        ads[k].obs[l] = k[i]
ad_merge = sc.concat(ads)

# preprocess
preprocessing.preprocess(
    ad_merge,
    normalize=is_raw,
    log1p=is_raw,
    scale="none",
    n_comps=n_comps,
    metric=metric,
    min_dist=min_dist,
    n_neighbors=n_neighbors,
    pca=True,
    umap=True,
    save_raw=False,
    min_counts=min_counts,
    min_genes=min_features,
    max_counts=max_counts,
    max_genes=max_features,
    min_cells=min_cells,
)

# save
df_umap = pd.DataFrame(ad_merge.obsm["X_umap"], index=ad_merge.obs_names, columns=["UMAP1", "UMAP2"])
df_umap[xenium_levels] = ad_merge.obs[xenium_levels]

df_umap.to_parquet(out_file)

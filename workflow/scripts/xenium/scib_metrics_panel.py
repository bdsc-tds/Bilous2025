import dask

dask.config.set({"dataframe.query-planning": False})

import argparse
import pandas as pd
import scanpy as sc
import scipy
from pathlib import Path
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from scib_metrics.benchmark._core import _LABELS


# Set up argument parser
parser = argparse.ArgumentParser(description="Embed panel of Xenium donors.")
parser.add_argument("--panel", type=Path, help="Path to the panel file.")
parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell type annotation dir.")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--normalisation", type=str, help="Normalisation method")
parser.add_argument("--layer", type=str, help="Name of saved layer of the seurat object, data or scale_data")
parser.add_argument("--reference", type=str, help="annotation reference")
parser.add_argument("--method", type=str, help="annotation method")
parser.add_argument("--level", type=str, help="annotation level")
parser.add_argument("--n_comps", type=int, help="Number of components.")
parser.add_argument("--max_n_cells", type=int, help="Max number of cells to use.")

args = parser.parse_args()

# Access the arguments
panel = args.panel
cell_type_annotation_dir = args.cell_type_annotation_dir
out_file = args.out_file
normalisation = args.normalisation
layer = args.layer
reference = args.reference
method = args.method
level = args.level
n_comps = args.n_comps
max_n_cells = args.max_n_cells

# variables
segmentation = panel.parents[1].stem
condition = panel.parents[0].stem
OBSM_KEY = "X_pca"
CT_KEY = (reference, method, level)
BATCH_KEY = "batch_key"
normalisation_annot = "lognorm"  # fix this for now, even for sctransfrom

# read xenium samples
ads = {}
for donor in (donors := panel.iterdir()):
    for sample in (samples := donor.iterdir()):
        k = (
            segmentation,
            condition,
            panel.stem,
            donor.stem,
            sample.stem,
        )
        name = "/".join(k)

        sample_counts_path = sample / f"{normalisation}/normalised_counts/{layer}.parquet"
        sample_idx_path = sample / f"{normalisation}/normalised_counts/cells.parquet"

        ads[k] = sc.AnnData(pd.read_parquet(sample_counts_path))
        if layer != "scale_data":  # no need to sparsify scale_data which is dense
            ads[k].X = scipy.sparse.csr_matrix(ads[k].X)
        ads[k].obs_names = pd.read_parquet(sample_idx_path).iloc[:, 0]

        # read cell type annotation
        sample_annotation_dir = cell_type_annotation_dir / f"{name}/{normalisation_annot}/reference_based"
        # for reference in (references := sample_annotation_dir.iterdir()):
        #     for method in (methods := reference.iterdir()):
        #         for level in (levels := method.iterdir()):

        annot_file = sample_annotation_dir / f"{reference}/{method}/{level}/single_cell/labels.parquet"
        ads[k].obs[(reference, method, level)] = pd.read_parquet(annot_file).set_index("cell_id").iloc[:, 0]


# concatenate
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
for k in ads.keys():
    for i, lvl in enumerate(xenium_levels):
        ads[k].obs[lvl] = k[i]
ad_merge = sc.concat(ads)
ad_merge.obs[BATCH_KEY] = ad_merge.obs[xenium_levels].agg("_".join, axis=1)

# drop NaN annotations
ad_merge = ad_merge[ad_merge.obs.notna().all(1)].copy()
# CT_KEYS = [c for c in ad_merge.obs.columns if c not in xenium_levels]

# subsample to reasonable size
if len(ad_merge) > max_n_cells:
    sc.pp.subsample(ad_merge, n_obs=max_n_cells)

# compute pca
sc.tl.pca(ad_merge, n_comps=n_comps)

# set up metrics
batchcor = BatchCorrection(
    silhouette_batch=True,
    ilisi_knn=True,
    kbet_per_label=False,
    graph_connectivity=False,
    pcr_comparison=False,
)

biocons = BioConservation(
    isolated_labels=False,
    nmi_ari_cluster_labels_leiden=True,
    nmi_ari_cluster_labels_kmeans=False,
    silhouette_label=True,
    clisi_knn=True,
)

# benchmark all cell type keys
df_results = pd.DataFrame()
# for i, CT_KEY in enumerate(CT_KEYS):
#     if i == 0:
bm = Benchmarker(
    ad_merge,
    batch_key=BATCH_KEY,
    label_key=CT_KEY,
    embedding_obsm_keys=[OBSM_KEY],
    pre_integrated_embedding_obsm_key=OBSM_KEY,
    bio_conservation_metrics=biocons,
    batch_correction_metrics=batchcor,
    n_jobs=-1,
)
bm.benchmark()
# else:
# to avoid recomputing kNN graph
# bm._emb_adatas[OBSM_KEY].obs[_LABELS] = ad_merge.obs[CT_KEY].values
# bm.benchmark()

# df_results[CT_KEY] = bm.get_results(min_max_scale=False).iloc[0]

# save results
bm.get_results(min_max_scale=False).iloc[[0]].to_parquet(out_file)

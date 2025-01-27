from pathlib import Path
import sys
import pandas as pd
import scanpy as sc

sys.path.append("workflow/scripts/")
import readwrite
import integration

panel_path = Path(sys.argv[1])
out_file = sys.argv[2]
n_comps = sys.argv[3]
n_neighbors = sys.argv[4]
metric = sys.argv[5]
min_dist = sys.argv[6]


# read xenium samples
xenium_paths = {}
for sample in (samples := panel_path.iterdir()):
    for replicate in (replicates := sample.iterdir()):
        k = (sample.stem, replicate.stem)
        replicate_path = replicate / "normalised_results/outs"
        name = "/".join(k)

        xenium_paths[k] = replicate_path

ads = readwrite.read_xenium_samples(
    xenium_paths, anndata_only=True, transcripts=False, sample_name_as_key=False
)


ad_merge = sc.concat(ads, label="dataset_merge_id")
integration.preprocess(
    ad_merge,
    normalize=True,
    log1p=True,
    scale="none",
    n_comps=n_comps,
    metric=metric,
    min_dist=min_dist,
    n_neighbors=n_neighbors,
    pca=True,
    umap=True,
    save_raw=False,
    filter_empty=True,
)

df_umap = pd.DataFrame(
    ad_merge.obsm["X_umap"], index=ad_merge.obs_names, columns=["UMAP1", "UMAP2"]
)
df_umap["dataset_merge_id"] = ad_merge.obs["dataset_merge_id"]
df_umap.to_parquet(out_file)

import pandas as pd
import numpy as np
import scanpy as sc
import sklearn
import spatialdata_io
import sys

# parameters
path = sys.argv[1]
out_file = sys.argv[2]
max_sample_size = sys.argv[3]
seed = 0

# read counts
adata = spatialdata_io.xenium(
    path,
    cells_as_circles=False,
    cells_boundaries=False,
    nucleus_boundaries=False,
    cells_labels=False,
    nucleus_labels=False,
    transcripts=False,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
)["table"]


# read annotations
annotations = {}
for reference in (
    references := path.parents[1] / "cell_type_annotation/reference_based"
).iterdir():
    for method in (methods := reference.iterdir()):
        for level in (levels := method.iterdir()):
            annotations[reference.stem, method.stem, level.stem] = pd.read_csv(
                level / "single_cell/labels.csv", index_col=0
            ).iloc[:, 0]

annotations = pd.DataFrame(annotations).dropna()
annotations.columns = ["_".join(col) for col in annotations.columns]

# add annotations to adata
adata = adata[annotations.index]
adata.obs = adata.obs.join(annotations)

# preprocess data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)

# precompute distances on subsample
sample_size = min(max_sample_size, len(adata))
indices = np.random.default_rng(seed).permutation(len(adata))[:sample_size]
D = sklearn.metrics.pairwise_distances(adata.obsm["X_pca"][indices])

# compute silhouettes
CT_KEYS = annotations.columns
silhouettes = pd.Series(index=CT_KEYS)
for CT_KEY in CT_KEYS:
    print(CT_KEY)
    labels = adata.obs[CT_KEY].iloc[indices].values
    silhouettes.loc[CT_KEY] = sklearn.metrics.silhouette_score(
        D, labels, random_state=seed, metric="precomputed"
    )

# save
silhouettes.to_frame().to_parquet(out_file)

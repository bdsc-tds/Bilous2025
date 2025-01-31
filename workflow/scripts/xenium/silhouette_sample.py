from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import sklearn
import spatialdata_io
import sys

# parameters
path = Path(sys.argv[1])
out_file = sys.argv[2]
max_donor_size = int(sys.argv[3])
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
adata.obs_names = adata.obs["cell_id"]

# read annotations
annotations = {}
for reference in (
    references := path.parents[1] / "cell_type_annotation/reference_based"
).iterdir():
    for method in (methods := reference.iterdir()):
        for level in (levels := method.iterdir()):
            if (level / "single_cell/labels.parquet").exists():
                annotations[reference.stem, method.stem, level.stem] = (
                    pd.read_parquet(level / "single_cell/labels.parquet")
                    .set_index("cell_id")
                    .iloc[:, 0]
                )


annotations = pd.DataFrame(annotations).dropna()
annotations.columns = [col for col in annotations.columns]

# add annotations to adata
adata = adata[annotations.index]
adata.obs = adata.obs.join(annotations)

# preprocess data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)

# precompute distances on subdonor
donor_size = min(max_donor_size, len(adata))
indices = np.random.default_rng(seed).permutation(len(adata))[:donor_size]
D = sklearn.metrics.pairwise_distances(adata.obsm["X_pca"][indices])

# compute silhouettes
CT_KEYS = annotations.columns
silhouettes = {}
for CT_KEY in CT_KEYS:
    print(CT_KEY)
    labels = adata.obs[CT_KEY].iloc[indices].values
    silhouettes[CT_KEY] = sklearn.metrics.silhouette_score(
        D, labels, random_state=seed, metric="precomputed"
    )

# save
silhouettes = pd.Series(silhouettes).reset_index()
silhouettes.columns = ["reference", "method", "level", "silhouette"]
silhouettes.to_parquet(out_file)

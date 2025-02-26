from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import sklearn
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Get silhouettes of a Xenium sample for all annotations.")
parser.add_argument("--sample_pca", type=str, help="Path to the sample PCA file.")
parser.add_argument("--sample_idx", type=str, help="Path to the sample index file.")
parser.add_argument("--sample_annotation_dir", type=Path, help="Path to the sample cell type annotation directory.")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--max_sample_size", type=int, help="Max number of samples to compute distance matrix")
parser.add_argument("--metric", type=str, help="Metric to compute distance matrix")

args = parser.parse_args()

# Access the arguments
sample_pca = args.sample_pca
sample_idx = args.sample_idx
sample_annotation_dir = args.sample_annotation_dir
out_file = args.out_file
max_sample_size = args.max_sample_size
metric = args.metric
seed = 0

# read counts
adata = sc.AnnData(pd.read_parquet(sample_pca))
adata.obs_names = pd.read_parquet(sample_idx).iloc[:, 0]

# read annotations
annotations = {}
for reference in (references := sample_annotation_dir.iterdir()):
    for method in (methods := reference.iterdir()):
        for level in (levels := method.iterdir()):
            if (level / "single_cell/labels.parquet").exists():
                annotations[reference.stem, method.stem, level.stem] = (
                    pd.read_parquet(level / "single_cell/labels.parquet").set_index("cell_id").iloc[:, 0]
                )


annotations = pd.DataFrame(annotations).dropna()
annotations.columns = [col for col in annotations.columns]

# add annotations to adata
adata = adata[annotations.index]
adata.obs = adata.obs.join(annotations)

# precompute distances on subsample
sample_size = min(max_sample_size, len(adata))
indices = np.random.default_rng(seed).permutation(len(adata))[:sample_size]
D = sklearn.metrics.pairwise_distances(adata.X[indices], metric=metric)

# compute silhouettes
CT_KEYS = annotations.columns
silhouettes = {}
for CT_KEY in CT_KEYS:
    print(CT_KEY)
    labels = adata.obs[CT_KEY].iloc[indices].values
    silhouettes[CT_KEY] = sklearn.metrics.silhouette_score(D, labels, random_state=seed, metric="precomputed")

# save
silhouettes = pd.Series(silhouettes).reset_index()
silhouettes.columns = ["reference", "method", "level", "silhouette"]
silhouettes.to_parquet(out_file)

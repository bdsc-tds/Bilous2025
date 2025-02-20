from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score

import dask

dask.config.set({"dataframe.query-planning": False})

import numpy as np
import pandas as pd
import squidpy as sq
import argparse
import os
import sys

sys.path.append("workflow/scripts/")
import preprocessing
import readwrite


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Run RESOLVI on a Xenium sample.")
    parser.add_argument("--path", type=str, help="Path to the xenium sample file.")
    parser.add_argument("--out_dir_resolvi_model", type=str, help="output directory with RESOLVI model weights")
    parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
    parser.add_argument(
        "--max_iter",
        type=int,
        default=50,
        help="Maximum number of iterations to train the model.",
    )
    parser.add_argument("--cell_type_labels", type=str, help="optional cell_type_labels for semi-supervised mode")
    parser.add_argument("--cti", type=int, help="cell type i")
    parser.add_argument("--ctj", type=int, help="cell type j")

    ret = parser.parse_args()
    if not os.path.isdir(ret.path):
        raise RuntimeError(f"Error! Input directory does not exist: {ret.path}")
    if ret.max_epochs <= 0:
        ret.max_epochs = 50

    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    # read counts
    adata = readwrite.read_xenium_sample(
        args.path,
        cells_as_circles=False,
        cells_boundaries=False,
        cells_boundaries_layers=False,
        nucleus_boundaries=False,
        cells_labels=False,
        nucleus_labels=False,
        transcripts=False,
        morphology_mip=False,
        morphology_focus=False,
        aligned_images=False,
        anndata=True,
    )

    # need to round proseg expected counts for resolVI to run
    # no need for if statement, doesn't change anything to other segmentation methods
    adata.X.data = adata.X.data.astype(np.float32).round()
    adata.obs_names = adata.obs_names.astype(str)

    labels_key = "labels_key"
    adata.obs[labels_key] = pd.read_parquet(args.cell_type_labels).iloc[:, 0]

    # preprocess (QC filters only)
    # resolvi requires at least 5 counts in each cell
    preprocessing.preprocess(
        adata,
        normalize=False,
        log1p=False,
        scale="none",
        pca=False,
        umap=False,
        save_raw=False,
        min_counts=args.min_counts,
        min_genes=args.min_features,
        max_counts=args.max_counts,
        max_genes=args.max_features,
        min_cells=args.min_cells,
        backend="cpu",
    )

    # Filter for neutrophils
    adata_cti = adata[adata.obs[labels_key] == args.cti]
    

    # Assume 'has_malignant_neighbor' is an annotation in adata.obs
    X = adata_cti.X  # Gene expression matrix
    y = adata_cti.obs["has_malignant_neighbor"]

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train logistic regression model
    model = LogisticRegression(max_iter=1000)
    model.fit(X_train, y_train)

    # Predict and evaluate
    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)

    # Inspect top informative genes
    importance = model.coef_[0]
    top_genes = neutrophils.var.index[importance.argsort()[-7:][::-1]]

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

import dask

dask.config.set({"dataframe.query-planning": False})

import scanpy as sc
import scipy
import numpy as np
import pandas as pd
import os
import sys
import argparse
import json
from pathlib import Path


sys.path.append("../../../workflow/scripts/")
import _utils
import readwrite


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Run RESOLVI on a Xenium sample.")
    parser.add_argument("--path", type=str, help="Path to the xenium sample file.")
    parser.add_argument("--out_dir_resolvi_model", type=str, help="output directory with RESOLVI model weights")
    parser.add_argument("--cell_type_labels", type=str, help="optional cell_type_labels for semi-supervised mode")
    parser.add_argument("--cti", type=int, help="cell type i")
    parser.add_argument("--ctj", type=int, help="cell type j")
    parser.add_argument("--sample_xeniumdir", type=str, help="")
    parser.add_argument("--sample_counts", type=str, help="")
    parser.add_argument("--sample_idx", type=str, help="")
    parser.add_argument("--cell_type_labels", type=str, help="")
    parser.add_argument("--out_file_permutation_summary", type=str, help="")
    parser.add_argument("--out_file_importances", type=str, help="")
    parser.add_argument("--out_file_importances_markers_rank_significance", type=str, help="")
    parser.add_argument("--n_neighbors", type=str, help="")
    parser.add_argument("--n_permutations", type=str, help="")
    parser.add_argument("--n_repeats", type=str, help="")
    parser.add_argument("--top_n", type=str, help="")
    parser.add_argument("--cti", type=str, help="")
    parser.add_argument("--ctj", type=str, help="")
    parser.add_argument("--scoring", type=str, help="")
    parser.add_argument("--markers", type=str, help="")

    ret = parser.parse_args()
    if not os.path.isdir(ret.path):
        raise RuntimeError(f"Error! Input directory does not exist: {ret.path}")

    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    # read raw data to get spatial coordinates
    adata = readwrite.read_xenium_sample(
        args.sample_dir,
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

    # read normalised data
    X_normalised = pd.read_parquet(args.sample_counts)
    X_normalised.index = pd.read_parquet(args.sample_idx).iloc[:, 0]
    adata = adata[X_normalised.index]
    adata.X = X_normalised

    # read labels
    label_key = "label_key"
    adata.obs[label_key] = pd.read_parquet(args.cell_type_labels).set_index("cell_id").iloc[:, 0]
    adata.obs[label_key] = adata.obs[label_key].replace(r" of .+", "", regex=True)

    # define target (cell type j presence in kNN)
    knnlabels = general_utils.get_knn_labels(adata, n_neighbors=args.n_neighbors, label_key=label_key, obsm="spatial")
    adata.obs[f"has_{args.ctj}_neighbor"] = knnlabels[args.ctj] > 0

    # read markers
    if args.markers == "diffexpr":
        sc.tl.rank_genes_groups(adata, groupby=label_key, groups=[args.ctj], reference="rest", method="wilcoxon")
        ctj_marker_genes = sc.get.rank_genes_groups_df(adata, group=args.ctj)[: args.top_n]
    else:
        df_markers = pd.read_json(args.markers)["canonical"].explode().reset_index()
        df_markers.columns = ["cell_type", "gene"]
        ctj_marker_genes = df_markers[df_markers["cell_type"] == args.ctj]["gene"].tolist()

    assert len(ctj_marker_genes), f"{args.ctj} not found in marker list or in DE list"

    # Filter for cti
    if args.cti is None:
        adata_cti = adata
    else:
        adata_cti = adata[adata.obs[label_key] == args.cti]

    # Split data
    X = adata_cti.X
    y = adata_cti.obs[f"has_{args.ctj}_neighbor"]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=0)

    # Init logistic regression model
    model = LogisticRegression()

    # Empirical p-value calculation using permutation test
    score, permutation_scores, p_value = permutation_test_score(
        model, X_train, y_train, scoring=args.scoring, n_permutations=args.n_permutations, n_jobs=-1
    )

    permutation_summary = pd.DataFrame(
        [[score, permutation_scores.mean(), permutation_scores.std(), p_value]],
        columns=[
            f"{args.scoring}_score",
            f"permutation_mean_{args.scoring}_score",
            f"permutation_std_{args.scoring}_score",
            "p_value",
        ],
    )
    permutation_summary["effect_size"] = (
        permutation_summary[f"{args.scoring}_score"] - permutation_summary[f"permutation_mean_{args.scoring}_score"]
    ) / permutation_summary["permutation_std_f1_score"]

    # Feature importances from permutations
    model.fit(X_train, y_train)
    importances = permutation_importance(
        model,
        pd.DataFrame.sparse.from_spmatrix(X_test),
        y_test,
        scoring=args.scoring,
        n_repeats=args.n_repeats,
        n_jobs=-1,
    )

    # Feature importances from model coefs
    # cv_results = cross_validate(model,X,y,return_estimator=True, scoring=scoring, n_jobs=-1)
    # importances = np.std(X, axis=0) * np.vstack([m.coef_[0] for m in cv_results['estimator']])

    # coef pvalues from formula
    # importances['pvalues'] = general_utils.logit_pvalue(model,X_train.toarray())[1:]

    # convert importances to df
    importances.pop("importances")
    importances = pd.DataFrame(importances, index=adata_cti.var_names).sort_values("importances_mean", ascending=False)

    # ctj marker rank significance from prerank
    markers_rank_significance = gseapy.prerank(
        rnk=importances["importances_mean"],
        gene_sets=[{"markers": ctj_marker_genes}],
        min_size=0,
    ).results

    # ctj marker rank significance from hypergeometric test
    N = len(importances)  # Total genes in ranked list
    K = len(ctj_marker_genes)  # Genes in the pathway/set of interest
    n = args.top_n  # Top ranked genes
    x = np.isin(importances.index[:n], ctj_marker_genes).sum()  # Overlapping genes in top n
    markers_rank_significance["hypergeometric_pvalue"] = scipy.stats.hypergeom.sf(x - 1, N, K, n)

    # Save
    permutation_summary.to_parquet(args.out_file_permutation_summary)
    importances.to_frame().to_parquet(args.out_file_importances)
    with open(args.out_file_importances_markers_rank_significance, "w") as f:
        json.dump(markers_rank_significance, f)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

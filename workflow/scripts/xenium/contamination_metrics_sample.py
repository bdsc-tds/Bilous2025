import dask

dask.config.set({"dataframe.query-planning": False})

import scanpy as sc
import pandas as pd
import argparse
import os
import numpy as np
import sys

sys.path.append("workflow/scripts/")
import _utils
import readwrite


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Run RESOLVI on a Xenium sample.")
    parser.add_argument("--sample_dir", type=str, help="")
    parser.add_argument("--sample_normalised_counts", type=str, help="")
    parser.add_argument("--sample_idx", type=str, help="")
    parser.add_argument("--sample_annotation", type=str, help="")
    parser.add_argument(
        "--out_file_df_permutations_logreg", type=str, help="path to the logreg permutation output file"
    )
    parser.add_argument("--out_file_df_importances_logreg", type=str, help="path to the logref importances output file")
    parser.add_argument("--out_file_df_diffexpr", type=str, help="path to the differential expression output file")
    parser.add_argument(
        "--out_file_df_markers_rank_significance_logreg",
        type=str,
        help="path to the logreg rank significance output file",
    )
    parser.add_argument(
        "--out_file_df_markers_rank_significance_diffexpr",
        type=str,
        help="path the differential expression rank significance output file",
    )
    # parser.add_argument("--out_dir_liana_lrdata", type=str, help="path to the liana lrdata results dir output")
    parser.add_argument("--n_neighbors", type=int, help="n° of neighbors to use to define the spatial graph")
    parser.add_argument("--n_permutations", type=int, help="n° of permutations for logreg random prediction baseline")
    parser.add_argument("--n_repeats", type=int, help="n° of repeats for logreg feature importances by permutations")
    parser.add_argument("--top_n", type=int, help="n° of top genes to evaluate for hypergeometric test")
    parser.add_argument(
        "--top_n_lr", type=int, help="n° of top ligand-receptor pairs to evaluate for liana hypergeometric test"
    )
    # parser.add_argument("--cti", type=str, help="cell type i (potentially corrupted by cell type j)")
    # parser.add_argument("--ctj", type=str, help="cell type j (potentially corrupting cell type i)")
    parser.add_argument("--scoring", type=str, help="sklearn scoring metric to use for logreg")
    parser.add_argument(
        "--markers",
        type=str,
        help="'diffexpr' to use empirical markers by diffexpr test, or a path to a df with cell_type and marker genes columns",
    )

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

    ####
    #### READ DATA
    ####
    # read raw counts and spatial coordinates
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
    if "proseg_expected" in args.sample_normalised_counts.as_posix():
        if "raw_results" not in args.sample_dir.as_posix():
            raise ValueError("raw results folder needed for proseg expected counts")
        adata.obs_names = "proseg-" + adata.obs_names.astype(str)

    # read normalised data, filter cells
    X_normalised = pd.read_parquet(args.sample_normalised_counts)
    X_normalised.index = pd.read_parquet(args.sample_idx).iloc[:, 0]
    X_normalised.columns = X_normalised.columns.str.replace(".", "-")  # undo seurat renaming
    adata = adata[X_normalised.index, X_normalised.columns]
    adata.layers["X_normalised"] = X_normalised

    # read labels
    label_key = "label_key"
    adata.obs[label_key] = pd.read_parquet(args.sample_annotation).set_index("cell_id").iloc[:, 0]
    adata.obs[label_key] = adata.obs[label_key].replace(r" of .+", "", regex=True)
    adata = adata[adata.obs[label_key].notna()]  # remove NaN annotation

    # read markers if needed
    if args.markers != "diffexpr":
        df_markers = pd.read_json(args.markers)["canonical"].explode().reset_index()
        df_markers.columns = ["cell_type", "gene"]

    # log-normalize before DE
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # define target (cell type j presence in kNN)
    obsm = "spatial"
    knnlabels, knndis, knnidx, knn_graph = _utils.get_knn_labels(
        adata, n_neighbors=args.n_neighbors, label_key=label_key, obsm=obsm, return_sparse_neighbors=True
    )

    adata.obsp[f"{args.obsm}_connectivities"] = knn_graph
    u_cell_types = adata.obs[label_key].unique()

    # iterate over cell types permutations (cell type i with cell type j presence in kNN)
    df_permutations_logreg = {}
    df_importances_logreg = {}
    df_diffexpr = {}
    df_markers_rank_significance_logreg = {}
    df_markers_rank_significance_diffexpr = {}
    df_markers_rank_significance_lrdata = {}

    for ctj in u_cell_types:
        if (adata.obs[label_key] == ctj).sum() < 30:
            print(f"Not enough cells from class {ctj}")
            continue

        # get markers
        if args.markers == "diffexpr":
            sc.tl.rank_genes_groups(adata, groupby=label_key, groups=[ctj], reference="rest", method="wilcoxon")
            ctj_marker_genes = sc.get.rank_genes_groups_df(adata, group=ctj)["names"][: args.top_n].tolist()
        else:
            ctj_marker_genes = df_markers[df_markers["cell_type"] == ctj]["gene"].tolist()
            ctj_marker_genes = [g for g in ctj_marker_genes if g in adata.var_names]

            assert len(ctj_marker_genes), f"no markers found for {ctj}"

        for cti in u_cell_types:
            if cti == ctj:
                continue
            print(cti, ctj)

            adata.obs[f"has_{ctj}_neighbor"] = knnlabels[ctj] > 0

            # read markers
            if args.markers == "diffexpr":
                sc.tl.rank_genes_groups(adata, groupby=label_key, groups=[ctj], reference="rest", method="wilcoxon")
                ctj_marker_genes = sc.get.rank_genes_groups_df(adata, group=ctj)[: args.top_n]
            else:
                ctj_marker_genes = df_markers[df_markers["cell_type"] == ctj]["gene"].tolist()
                ctj_marker_genes = [g for g in ctj_marker_genes if g in adata.var_names]

            assert len(ctj_marker_genes), f"{ctj} not found in marker list or in DE list"

            # Filter for cti
            adata_cti = adata[adata.obs[label_key] == cti]

            if (adata_cti.obs[f"has_{ctj}_neighbor"]).sum() < 30 or (~adata_cti.obs[f"has_{ctj}_neighbor"]).sum() < 30:
                print(f"Not enough cells from each class to test {cti} with {ctj} neighbors")
                continue

            ####
            #### LOGISTIC REGRESSION TEST: predict ctj in kNN based on cti expression
            ####

            # train logreg model
            df_permutations_logreg[cti, ctj], df_importances_logreg[cti, ctj] = _utils.logreg(
                X=adata_cti.layers["X_normalised"],
                y=adata_cti.obs[f"has_{ctj}_neighbor"],
                feature_names=adata.var_names,
                scoring=args.scoring,
                test_size=0.2,
                n_permutations=args.n_permutations,
                n_repeats=args.n_repeats,
                random_state=0,
                max_iter=1000,
                importance_mode="coef",
            )

            # get significance from gsea and hypergeometric test
            df_markers_rank_significance_logreg[cti, ctj] = _utils.get_marker_rank_significance(
                rnk=df_importances_logreg[cti, ctj]["importances"], gene_set=ctj_marker_genes, top_n=args.top_n
            )

            ###
            ### DIFF EXPR TEST: check DE genes between cti with ctj neighbor or not
            ###
            adata_cti.obs[f"has_{ctj}_neighbor_str"] = adata_cti.obs[f"has_{ctj}_neighbor"].astype(str)
            sc.tl.rank_genes_groups(
                adata_cti, groupby=f"has_{ctj}_neighbor_str", groups=["True"], reference="False", method="wilcoxon"
            )
            df_diffexpr[cti, ctj] = sc.get.rank_genes_groups_df(adata_cti, group="True")
            df_diffexpr[cti, ctj]["rank_score"] = (
                -np.log10(df_diffexpr[cti, ctj]["pvals_adj"]) * df_diffexpr[cti, ctj]["logfoldchanges"]
            )
            df_diffexpr[cti, ctj] = df_diffexpr[cti, ctj].sort_values("rank_score", ascending=False)

            # get significance from gsea and hypergeometric test
            df_markers_rank_significance_diffexpr[cti, ctj] = _utils.get_marker_rank_significance(
                rnk=df_diffexpr[cti, ctj].set_index("names")["rank_score"],
                gene_set=ctj_marker_genes,
                top_n=args.top_n,
            )

    ###
    ### CONCAT AND SAVE OUTPUTS
    ###
    # logreg
    pd.concat(df_permutations_logreg).to_parquet(args.out_file_df_permutations_logreg)
    pd.concat(df_importances_logreg).to_parquet(args.out_file_df_importances_logreg)
    pd.concat(df_markers_rank_significance_logreg).to_parquet(args.out_file_df_markers_rank_significance_logreg)

    # diffexpr
    pd.concat(df_diffexpr).to_parquet(args.out_file_df_diffexpr)
    pd.concat(df_markers_rank_significance_diffexpr).to_parquet(args.out_file_df_markers_rank_significance_diffexpr)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

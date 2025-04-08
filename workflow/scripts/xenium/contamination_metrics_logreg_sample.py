import dask

dask.config.set({"dataframe.query-planning": False})

import numpy as np
import scanpy as sc
import pandas as pd
import argparse
import os
import sys

sys.path.append("workflow/scripts/")
import _utils
import readwrite
import preprocessing


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Run RESOLVI on a Xenium sample.")
    parser.add_argument("--sample_corrected_counts_path", type=str, help="")
    parser.add_argument("--sample_dir", type=str, help="")
    parser.add_argument("--sample_normalised_counts", type=str, help="")
    parser.add_argument("--sample_idx", type=str, help="")
    parser.add_argument("--sample_annotation", type=str, help="")
    parser.add_argument(
        "--out_file_df_permutations_logreg", type=str, help="path to the logreg permutation output file"
    )
    parser.add_argument("--out_file_df_importances_logreg", type=str, help="path to the logref importances output file")
    parser.add_argument(
        "--out_file_df_markers_rank_significance_logreg",
        type=str,
        help="path to the logreg rank significance output file",
    )
    parser.add_argument("--radius", type=int, help="n° of neighbors to use to define the spatial graph")
    parser.add_argument("--cv_mode", type=str, help="cv_mode for logreg")
    parser.add_argument("--n_splits", type=int, help="n° of splits for logreg random prediction baseline")
    parser.add_argument("--n_permutations", type=int, help="n° of permutations for logreg random prediction baseline")
    parser.add_argument("--n_repeats", type=int, help="n° of repeats for logreg feature importances by permutations")
    parser.add_argument("--top_n", type=int, help="n° of top genes to evaluate for hypergeometric test")
    parser.add_argument("--scoring", type=str, help="sklearn scoring metric to use for logreg")
    parser.add_argument(
        "--markers",
        type=str,
        help="'diffexpr' to use empirical markers by diffexpr test, or a path to a df with cell_type and marker genes columns",
    )
    parser.add_argument(
        "--precomputed_ctj_markers",
        type=str,
        help="path to ctj markers parquet obtained from running the function on raw counts",
    )
    parser.add_argument(
        "--precomputed_adata_obs",
        type=str,
        help="path to adata obs parquet obtained from running the function on raw counts",
    )
    parser.add_argument("--min_counts", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--min_features", type=int, help="QC parameter from pipeline config")
    parser.add_argument("--max_counts", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--max_features", type=float, help="QC parameter from pipeline config")
    parser.add_argument("--min_cells", type=int, help="QC parameter from pipeline config")
    parser.add_argument(
        "--max_n_cells",
        type=int,
        help="max_n_cells to use",
    )
    parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )
    ret = parser.parse_args()
    return ret


if __name__ == "__main__":
    args = parse_args()

    if args.l is not None:
        old_stdout = sys.stdout
        old_stderr = sys.stderr

        _log = open(args.l, "w", encoding="utf-8")

        sys.stdout = _log
        sys.stderr = _log

    list_n_markers = [10, 20, 30, 40, 50]
    obsm = "spatial"
    label_key = "label_key"
    index_diffexpr_metrics = [
        "Name",
        "Term",
        "ES",
        "NES",
        "NOM p-val",
        "FDR q-val",
        "FWER p-val",
        "Tag %",
        "Gene %",
        "Lead_genes",
        "hypergeometric_pvalue",
        f"n_hits_top_n={args.top_n}",
        "mean_zscore",
        "mean_zscore_pvalue",
    ]

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
    if "proseg_expected" in args.sample_normalised_counts:
        if "raw_results" not in args.sample_dir:
            raise ValueError("raw results folder needed for proseg expected counts")
        adata.obs_names = "proseg-" + adata.obs_names.astype(str)

    # read corrected counts
    if args.sample_corrected_counts_path is not None:
        adata_corrected_counts = sc.read_10x_h5(
            args.sample_corrected_counts_path,
        )

        adata_corrected_counts.obsm[obsm] = adata[adata_corrected_counts.obs_names].obsm[obsm]
        adata = adata_corrected_counts

    # read normalised data, filter cells
    # X_normalised = pd.read_parquet(args.sample_normalised_counts)
    # X_normalised.index = pd.read_parquet(args.sample_idx).iloc[:, 0]
    # X_normalised.columns = X_normalised.columns.str.replace(".", "-")  # undo seurat renaming
    # obs_found = [c for c in X_normalised.index if c in adata.obs_names]
    # var_found = [g for g in X_normalised.columns if g in adata.var_names]
    # adata = adata[obs_found, var_found]
    # adata.layers["X_normalised"] = X_normalised.loc[obs_found, var_found]

    # reapply QC to corrected counts data
    preprocessing.preprocess(
        adata,
        min_counts=args.min_counts,
        min_genes=args.min_features,
        max_counts=args.max_counts,
        max_genes=args.max_features,
        min_cells=args.min_cells,
        save_raw=False,
    )

    # read labels
    label_key = "label_key"
    adata.obs[label_key] = pd.read_parquet(args.sample_annotation).set_index("cell_id").iloc[:, 0]
    adata = adata[adata.obs[label_key].notna()]  # remove NaN annotation

    if "Level2.1" in args.sample_annotation:
        # for custom Level2.1, simplify subtypes
        adata.obs.loc[adata.obs[label_key].str.contains("malignant"), label_key] = "malignant cell"
        adata.obs.loc[adata.obs[label_key].str.contains("T cell"), label_key] = "T cell"

    # read markers if needed
    if args.markers != "diffexpr":
        if "Level2.1" in args.sample_annotation and args.markers == "common_markers":
            # custom mapping for Level2.1: simplify to Level1 to assess with known markers
            level_simplified = "Level1"

            df_markers = pd.read_csv("data/markers/Xenium_panels_common_markers.csv")[["cell_type", "gene_name"]]
            palette = pd.read_csv("data/xenium/metadata/col_palette_cell_types_combo.csv")

            cell_types_mapping = palette.set_index("Level2.1")[level_simplified].replace(r" of .+", "", regex=True)
            cell_types_mapping[cell_types_mapping.str.contains("malignant")] = "malignant cell"
            adata.obs[label_key] = adata.obs[label_key].replace(cell_types_mapping)

        else:
            df_markers = pd.read_csv(args.markers)[["cell_type", "gene_name"]]

        ct_not_found = adata.obs[label_key][~adata.obs[label_key].isin(df_markers["cell_type"])].unique()
        print(f"Could not find {ct_not_found} in markers file")
        adata = adata[adata.obs[label_key].isin(df_markers["cell_type"])]
    else:
        # get precomputed markers from raw data
        if args.precomputed_ctj_markers is not None:
            print("Loading precomputed markers")
            df_ctj_marker_genes_precomputed = pd.read_parquet(args.precomputed_ctj_markers)

    ####
    #### PREPARE DATA
    ####
    # log-normalize before DE
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # define target (cell type j presence in kNN)
    if args.precomputed_adata_obs is not None:
        print("Loading precomputed adata obs. Replacing loaded labels")
        adata.obs = pd.read_parquet(args.precomputed_adata_obs).loc[adata.obs_names]
    else:
        obsm = "spatial"
        knnlabels, knndis, knnidx, knn_graph = _utils.get_knn_labels(
            adata, radius=args.radius, label_key=label_key, obsm=obsm, return_sparse_neighbors=True
        )

        adata.obsp[f"{obsm}_connectivities"] = knn_graph

    ####
    #### SCORE CONTAMINATION
    ####
    # iterate over cell types permutations (cell type i with cell type j presence in kNN)
    u_cell_types = adata.obs[label_key].unique()

    df_permutations_logreg = {}
    df_importances_logreg = {}
    df_markers_rank_significance_logreg = {}

    for ctj in u_cell_types:
        if (adata.obs[label_key] == ctj).sum() < 30:
            print(f"Not enough cells from class {ctj}")
            continue

        # get markers
        if args.markers == "diffexpr":
            if args.precomputed_ctj_markers is not None:
                ctj_marker_genes_precomputed = df_ctj_marker_genes_precomputed[ctj]

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

            if args.precomputed_adata_obs is None:
                adata.obs[f"has_{ctj}_neighbor"] = knnlabels[ctj] > 0

            # Filter for cti
            adata_cti = adata[adata.obs[label_key] == cti]
            if len(adata_cti) > args.max_n_cells:
                rng = np.random.default_rng(0)
                print(f"Subsampling {cti} to {args.max_n_cells} cells")

                subsampled_true = adata_cti.obs_names[adata_cti.obs[f"has_{ctj}_neighbor"]]
                subsampled_true = rng.choice(
                    subsampled_true, min(len(subsampled_true), args.max_n_cells // 2), replace=False
                )
                subsampled_false = adata_cti.obs_names[~adata_cti.obs[f"has_{ctj}_neighbor"]]
                subsampled_false = rng.choice(
                    subsampled_false, min(len(subsampled_false), args.max_n_cells // 2), replace=False
                )

                subsampled_idx = np.hstack([subsampled_true, subsampled_false])
                adata_cti = adata_cti[subsampled_idx]

            if (adata_cti.obs[f"has_{ctj}_neighbor"]).sum() < 30 or (~adata_cti.obs[f"has_{ctj}_neighbor"]).sum() < 30:
                print(f"Not enough cells from each class to test {cti} with {ctj} neighbors")
                continue

            ####
            #### LOGISTIC REGRESSION TEST: predict ctj in kNN based on cti expression
            ####

            # train logreg model
            df_permutations_logreg[cti, ctj], df_importances_logreg[cti, ctj] = _utils.logreg(
                X=adata_cti.X,  # adata_cti.layers["X_normalised"],
                y=adata_cti.obs[f"has_{ctj}_neighbor"],
                feature_names=adata.var_names,
                scoring=args.scoring,
                test_size=0.2,
                n_splits=args.n_splits,
                n_permutations=args.n_permutations,
                n_repeats=args.n_repeats,
                random_state=0,
                max_iter=500,
                importance_mode="coef",
                class_weight="balanced",
                cv_mode=args.cv_mode,
                spatial_coords=adata_cti.obsm["spatial"],
            )

            # get significance from gsea and hypergeometric test
            rank_metric = "importances"
            df_markers_rank_significance_logreg[cti, ctj] = pd.DataFrame(index=index_diffexpr_metrics)
            dict_ctj_marker_genes = {"": ctj_marker_genes}

            if args.precomputed_ctj_markers is not None:
                # also compute scores for precomputed marker gene list
                dict_ctj_marker_genes["_precomputed"] = ctj_marker_genes_precomputed

            for n in list_n_markers:
                for k_, markers_ in dict_ctj_marker_genes.items():
                    markers_n_ = markers_[:n]

                    df_markers_rank_significance_logreg[cti, ctj][rank_metric + k_ + f"_{n=}"] = (
                        _utils.get_marker_rank_significance(
                            rnk=df_importances_logreg[cti, ctj][rank_metric].sort_values(ascending=False),
                            gene_set=markers_n_,
                            top_n=args.top_n,
                        ).iloc[0]
                    )

    ###
    ### CONCAT AND SAVE OUTPUTS
    ###

    # logreg
    pd.concat(df_permutations_logreg).to_parquet(args.out_file_df_permutations_logreg)
    pd.concat(df_importances_logreg).to_parquet(args.out_file_df_importances_logreg)
    pd.concat(df_markers_rank_significance_logreg).T.to_parquet(args.out_file_df_markers_rank_significance_logreg)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

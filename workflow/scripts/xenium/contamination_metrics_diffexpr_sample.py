import dask

dask.config.set({"dataframe.query-planning": False})

import json
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import os
import scipy
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
        "--out_file_df_ctj_marker_genes",
        type=str,
        help="path to the marker genes used for each contaminating cell type",
    )
    parser.add_argument("--out_file_df_diffexpr", type=str, help="path to the differential expression output file")
    parser.add_argument(
        "--out_file_df_markers_rank_significance_diffexpr",
        type=str,
        help="path to the differential expression rank significance output file",
    )
    parser.add_argument("--out_file_summary_stats", type=str, help="path to the summary stats output file")
    parser.add_argument("--out_file_adata_obs", type=str, help="path to the adata.obs output file")
    parser.add_argument("--radius", type=int, help="n° of neighbors to use to define the spatial graph")
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

    obsm = "spatial"
    label_key = "label_key"
    rank_metrics = ["logfoldchanges", "-log10pvals_x_logfoldchanges", "-log10pvals_x_sign_logfoldchanges"]
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

        adata_corrected_counts.obsm["spatial"] = adata[adata_corrected_counts.obs_names].obsm["spatial"]
        adata = adata_corrected_counts

    # read normalised data, filter cells
    # X_normalised = pd.read_parquet(args.sample_normalised_counts)
    # X_normalised.index = pd.read_parquet(args.sample_idx).iloc[:, 0]
    # X_normalised.columns = X_normalised.columns.str.replace(".", "-")  # undo seurat renaming
    # obs_found = [c for c in X_normalised.index if c in adata.obs_names]
    # var_found = [g for g in X_normalised.columns if g in adata.var_names]
    # adata = adata[obs_found,var_found]
    # adata = adata[X_normalised.index, X_normalised.columns]
    # adata.layers["X_normalised"] = X_normalised.loc[obs_found,var_found]

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

    adata.obs[label_key] = pd.read_parquet(args.sample_annotation).set_index("cell_id").iloc[:, 0]
    adata = adata[adata.obs[label_key].notna()]  # remove NaN annotation

    if "Level2.1" in args.sample_annotation:
        # for custom Level2.1, simplify malignant subtypes to malignant
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
            print(f"Loading precomputed {args.ctj} markers")
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
        adata.obs = pd.read_parquet(args.precomputed_adata_obs)
    else:
        knnlabels, knndis, knnidx, knn_graph = _utils.get_knn_labels(
            adata, radius=args.radius, label_key=label_key, obsm=obsm, return_sparse_neighbors=True
        )

        adata.obsp[f"{obsm}_connectivities"] = knn_graph

    ####
    #### SCORE CONTAMINATION
    ####
    # iterate over cell types permutations (cell type i with cell type j presence in kNN)
    u_cell_types = adata.obs[label_key].unique()

    df_diffexpr = {}
    df_markers_rank_significance_diffexpr = {}
    df_ctj_marker_genes = {}

    for ctj in u_cell_types:
        if (adata.obs[label_key] == ctj).sum() < 30:
            print(f"Not enough cells from class {ctj}")
            continue

        # get markers
        if args.markers == "diffexpr":
            if args.precomputed_ctj_markers is not None:
                print(f"Loading precomputed {ctj} markers")
                ctj_marker_genes_precomputed = pd.read_parquet(args.precomputed_ctj_markers)

            sc.tl.rank_genes_groups(adata, groupby=label_key, groups=[ctj], reference="rest", method="wilcoxon")
            ctj_marker_genes = sc.get.rank_genes_groups_df(adata, group=ctj)["names"][: args.top_n].tolist()
        else:
            ctj_marker_genes = df_markers[df_markers["cell_type"] == ctj]["gene"].tolist()
            ctj_marker_genes = [g for g in ctj_marker_genes if g in adata.var_names]

            assert len(ctj_marker_genes), f"no markers found for {ctj}"

        df_ctj_marker_genes[ctj] = ctj_marker_genes

        for cti in u_cell_types:
            if cti == ctj:
                continue
            print(cti, ctj)

            adata.obs[f"has_{ctj}_neighbor"] = knnlabels[ctj] > 0

            # Filter for cti
            adata_cti = adata[adata.obs[label_key] == cti]

            if (adata_cti.obs[f"has_{ctj}_neighbor"]).sum() < 30 or (~adata_cti.obs[f"has_{ctj}_neighbor"]).sum() < 30:
                print(f"Not enough cells from each class to test {cti} with {ctj} neighbors")
                continue

            ###
            ### DIFF EXPR TEST: check DE genes between cti with ctj neighbor or not
            ###

            adata_cti.obs[f"has_{ctj}_neighbor_str"] = adata_cti.obs[f"has_{ctj}_neighbor"].astype(str)
            sc.tl.rank_genes_groups(
                adata_cti, groupby=f"has_{ctj}_neighbor_str", groups=["True"], reference="False", method="wilcoxon"
            )
            df_diffexpr[cti, ctj] = sc.get.rank_genes_groups_df(adata_cti, group="True")

            # add ranking score by -log10pvals_x_logfoldchanges. Use offset to avoid -log10(pval) = inf
            pvals = df_diffexpr[cti, ctj]["pvals"]
            df_diffexpr[cti, ctj]["pvals_offset"] = pvals + pvals[pvals > 0].min() * 0.1
            df_diffexpr[cti, ctj]["-log10pvals_x_logfoldchanges"] = (
                -np.log10(df_diffexpr[cti, ctj]["pvals_offset"]) * df_diffexpr[cti, ctj]["logfoldchanges"]
            )

            # add ranking score by log10pvals_x_signlogFC
            df_diffexpr[cti, ctj]["-log10pvals_x_sign_logfoldchanges"] = -np.log10(
                df_diffexpr[cti, ctj]["pvals"]
            ) * np.sign(df_diffexpr[cti, ctj]["logfoldchanges"])

            # get significance from gsea and hypergeometric test
            df_markers_rank_significance_diffexpr[cti, ctj] = pd.DataFrame(index=index_diffexpr_metrics)
            dict_ctj_marker_genes = {"": ctj_marker_genes}

            if args.precomputed_ctj_markers is not None:
                # also compute scores for precomputed marker gene list
                dict_ctj_marker_genes["_precomputed"] = ctj_marker_genes_precomputed

            for k_, markers_ in dict_ctj_marker_genes.items():
                for rank_metric in rank_metrics:
                    df_markers_rank_significance_diffexpr[cti, ctj][rank_metric + k_] = (
                        _utils.get_marker_rank_significance(
                            rnk=df_diffexpr[cti, ctj].set_index("names")[rank_metric].sort_values(ascending=False),
                            gene_set=markers_,
                            top_n=args.top_n,
                        ).iloc[0]
                    )

                # add pvalue metrics
                mean_zscore = df_diffexpr[cti, ctj].set_index("names")["scores"].loc[markers_].mean()
                mean_zscore_pvalue = scipy.stats.norm.sf(np.abs(mean_zscore)) * 2  # Two-tailed p-value
                df_markers_rank_significance_diffexpr[cti, ctj]["mean_zscore" + k_] = pd.Series(
                    [mean_zscore, mean_zscore_pvalue], index=["mean_zscore", "mean_zscore_pvalue"]
                )

    # count number of True/False for each has_{ctj}_neighbor column
    cols = [f"has_{ctj}_neighbor" for ctj in u_cell_types if f"has_{ctj}_neighbor" in adata.obs.columns]
    df_has_neighbor_counts = (
        adata.obs.melt(id_vars=[label_key], value_vars=cols)
        .groupby([label_key, "variable"], observed=True)["value"]
        .value_counts()
        .reset_index(name="count")
    )

    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1

    ###
    ### CONCAT AND SAVE OUTPUTS
    ###
    # general stats
    summary_stats = {
        "n_cells": len(adata),
        "n_cells_by_type": adata.obs[label_key].value_counts().to_dict(),
        "mean_n_genes_by_type": adata.obs.groupby(label_key, observed=True)["n_genes"].mean().to_dict(),
        "median_n_genes_by_type": adata.obs.groupby(label_key, observed=True)["n_genes"].median().to_dict(),
        "mean_n_genes": adata.obs["n_genes"].mean(),
        "median_n_genes": adata.obs["n_genes"].median(),
        "df_has_neighbor_counts": df_has_neighbor_counts.to_dict(),  # Storing the DataFrame
    }

    with open(args.out_file_summary_stats, "w") as f:
        json.dump(summary_stats, f)

    # obs
    adata.obs.to_parquet(args.out_file_adata_obs)
    # markers
    df_ctj_marker_genes = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in df_ctj_marker_genes.items()]))
    df_ctj_marker_genes.to_parquet(args.out_file_df_ctj_marker_genes)
    # diffexpr
    pd.concat(df_diffexpr).to_parquet(args.out_file_df_diffexpr)
    pd.concat(df_markers_rank_significance_diffexpr).T.to_parquet(args.out_file_df_markers_rank_significance_diffexpr)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

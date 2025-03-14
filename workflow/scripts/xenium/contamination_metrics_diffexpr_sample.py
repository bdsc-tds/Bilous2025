import dask

dask.config.set({"dataframe.query-planning": False})

import scanpy as sc
import pandas as pd
import argparse
import os
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
    parser.add_argument("--n_neighbors", type=int, help="n° of neighbors to use to define the spatial graph")
    parser.add_argument("--top_n", type=int, help="n° of top genes to evaluate for hypergeometric test")
    parser.add_argument("--scoring", type=str, help="sklearn scoring metric to use for logreg")
    parser.add_argument(
        "--markers",
        type=str,
        help="'diffexpr' to use empirical markers by diffexpr test, or a path to a df with cell_type and marker genes columns",
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

    # read normalised data, filter cells
    X_normalised = pd.read_parquet(args.sample_normalised_counts)
    X_normalised.index = pd.read_parquet(args.sample_idx).iloc[:, 0]
    X_normalised.columns = X_normalised.columns.str.replace(".", "-")  # undo seurat renaming
    adata = adata[X_normalised.index, X_normalised.columns]
    adata.layers["X_normalised"] = X_normalised

    # read labels
    label_key = "label_key"
    adata.obs[label_key] = pd.read_parquet(args.sample_annotation).set_index("cell_id").iloc[:, 0]
    adata = adata[adata.obs[label_key].notna()]  # remove NaN annotation

    if "Level2.1" in args.sample_annotation:
        # for custom Level2.1, simplify malignant subtypes to malignant
        adata.obs.loc[adata.obs[label_key].str.contains("malignant"), label_key] = "malignant cell"

    # read markers if needed
    if args.markers != "diffexpr":
        if "Level2.1" in args.sample_annotation:
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

    # log-normalize before DE
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # define target (cell type j presence in kNN)
    obsm = "spatial"
    knnlabels, knndis, knnidx, knn_graph = _utils.get_knn_labels(
        adata, n_neighbors=args.n_neighbors, label_key=label_key, obsm=obsm, return_sparse_neighbors=True
    )

    adata.obsp[f"{obsm}_connectivities"] = knn_graph
    u_cell_types = adata.obs[label_key].unique()

    # iterate over cell types permutations (cell type i with cell type j presence in kNN)
    df_diffexpr = {}
    df_markers_rank_significance_diffexpr = {}
    df_markers_rank_significance_lrdata = {}
    df_ctj_marker_genes = {}

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
            df_diffexpr[cti, ctj] = sc.get.rank_genes_groups_df(adata_cti, group="True").sort_values("pvals_adj")

            # get significance from gsea and hypergeometric test
            df_markers_rank_significance_diffexpr[cti, ctj] = _utils.get_marker_rank_significance(
                rnk=df_diffexpr[cti, ctj].set_index("names")["logfoldchanges"],
                gene_set=ctj_marker_genes,
                top_n=args.top_n,
            )

    ###
    ### CONCAT AND SAVE OUTPUTS
    ###
    # markers
    df_ctj_marker_genes = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in df_ctj_marker_genes.items()]))
    df_ctj_marker_genes.to_parquet(args.out_file_df_ctj_marker_genes)
    # diffexpr
    pd.concat(df_diffexpr).to_parquet(args.out_file_df_diffexpr)
    pd.concat(df_markers_rank_significance_diffexpr).to_parquet(args.out_file_df_markers_rank_significance_diffexpr)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

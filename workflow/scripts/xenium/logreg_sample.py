import dask

dask.config.set({"dataframe.query-planning": False})

# import liana
import scanpy as sc
import pandas as pd
import sys
import argparse
import os
from pathlib import Path

sys.path.append("../../../workflow/scripts/")
import _utils
import readwrite


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Run RESOLVI on a Xenium sample.")
    parser.add_argument("--sample_dir", type=str, help="")
    parser.add_argument("--sample_counts", type=str, help="")
    parser.add_argument("--sample_idx", type=str, help="")
    parser.add_argument("--cell_type_labels", type=str, help="")
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
    parser.add_argument("--out_dir_liana_lrdata", type=str, help="path to the liana lrdata results dir output")
    parser.add_argument("--n_neighbors", type=int, help="n° of neighbors to use to define the spatial graph")
    parser.add_argument("--n_permutations", type=int, help="n° of permutations for logreg random prediction baseline")
    parser.add_argument("--n_repeats", type=int, help="n° of repeats for logreg feature importances by permutations")
    parser.add_argument("--top_n", type=int, help="n° of top genes to evaluate for hypergeometric test")
    parser.add_argument("--cti", type=str, help="cell type i (potentially corrupted by cell type j)")
    parser.add_argument("--ctj", type=str, help="cell type j (potentially corrupting cell type i)")
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

    # read normalised data, filter cells
    X_normalised = pd.read_parquet(args.sample_counts)
    X_normalised.index = pd.read_parquet(args.sample_idx).iloc[:, 0]
    adata = adata[X_normalised.index]
    adata.X = X_normalised

    # read labels
    label_key = "label_key"
    adata.obs[label_key] = pd.read_parquet(args.cell_type_labels).set_index("cell_id").iloc[:, 0]
    adata.obs[label_key] = adata.obs[label_key].replace(r" of .+", "", regex=True)

    # define target (cell type j presence in kNN)
    obsm = "spatial"
    knnlabels, knndis, knnidx, knn_graph = _utils.get_knn_labels(
        adata, n_neighbors=args.n_neighbors, label_key=label_key, obsm=obsm, return_sparse_neighbors=True
    )

    adata.obsp[f"{args.obsm}_connectivities"] = knn_graph
    adata.obs[f"count_{args.ctj}_neighbor"] = knnlabels[args.ctj]
    adata.obs[f"has_{args.ctj}_neighbor"] = knnlabels[args.ctj] > 0

    # read markers
    if args.markers == "diffexpr":
        sc.tl.rank_genes_groups(adata, groupby=label_key, groups=[args.ctj], reference="rest", method="wilcoxon")
        ctj_marker_genes = sc.get.rank_genes_groups_df(adata, group=args.ctj)[: args.top_n]
    else:
        df_markers = pd.read_json(args.markers)["canonical"].explode().reset_index()
        df_markers.columns = ["cell_type", "gene"]
        ctj_marker_genes = df_markers[df_markers["cell_type"] == args.ctj]["gene"].tolist()
        ctj_marker_genes = [g for g in ctj_marker_genes if g in adata.var_names]

    assert len(ctj_marker_genes), f"{args.ctj} not found in marker list or in DE list"

    # Filter for cti
    if args.cti is None:
        adata_cti = adata
    else:
        adata_cti = adata[adata.obs[label_key] == cti]

    ####
    #### LOGISTIC REGRESSION TEST: predict ctj in kNN based on cti expression
    ####

    # train logreg model
    df_permutations_logreg, df_importances_logreg = _utils.logreg(
        X=adata_cti.X,
        y=adata_cti.obs[f"has_{args.ctj}_neighbor"],
        feature_names=adata.var_names,
        scoring=args.scoring,
        test_size=0.2,
        n_permutations=args.n_permutations,
        n_repeats=args.n_repeats,
        random_state=0,
    )

    # get significance from gsea and hypergeometric test
    df_markers_rank_significance_logreg = _utils.get_marker_rank_significance(
        rnk=df_importances_logreg["importances_mean"], gene_set=ctj_marker_genes, top_n=args.top_n
    )

    ###
    ### DIFF EXPR TEST: check DE genes between cti with ctj neighbor or not
    ###
    idx_no_ctj_neighbor = adata_cti.obs[f"count_{args.ctj}_neighbor"] == 0
    if sum(idx_no_ctj_neighbor) < 30:  # arbitrary threshold to consider there's enough cells for DE
        raise ValueError("Not enough cells without ctj neighbors")

    adata_cti.obs[f"has_{args.ctj}_neighbor_str"] = adata_cti.obs[f"has_{args.ctj}_neighbor"].astype(str)
    sc.tl.rank_genes_groups(
        adata_cti, groupby=f"has_{args.ctj}_neighbor_str", groups=["True"], reference="False", method="wilcoxon"
    )
    df_diffexpr = sc.get.rank_genes_groups_df(adata_cti, group="True").sort_values("pvals_adj")

    # get significance from gsea and hypergeometric test
    df_markers_rank_significance_diffexpr = _utils.get_marker_rank_significance(
        rnk=df_diffexpr.set_index("names")["logfoldchanges"], gene_set=ctj_marker_genes, top_n=args.top_n
    )

    ###
    ### CELL-CELL COMMUNICATION TEST: check communication between cti with ctj neighbor
    ###
    # adata = adata[adata.obs[f'has_{ctj}_neighbor']]
    # lrdata = liana.mt.bivariate(
    #     adata,
    #     connectivity_key = f'{obsm}_connectivities',
    #     resource_name='consensus', # NOTE: uses HUMAN gene symbols!
    #     local_name='cosine', # Name of the function
    #     global_name='morans',
    #     n_perms=30, # Number of permutations to calculate a p-value
    #     mask_negatives=True, # Whether to mask LowLow/NegativeNegative interactions
    #     add_categories=True, # Whether to add local categories to the results
    #     nz_prop=0.0, # Minimum expr. proportion for ligands/receptors and their subunits
    #     use_raw=False,
    #     verbose=True
    #     )

    # # get significance from gsea and hypergeometric test
    # df_markers_rank_significance_diffexpr = _utils.get_marker_rank_significance(
    #     rnk=df_diffexpr.set_index('names')['logfoldchanges'],
    #     gene_set=ctj_marker_genes,
    #     top_n = top_n)

    ###
    ### SAVE OUTPUTS
    ###
    # logreg
    df_permutations_logreg.to_parquet(args.out_file_df_permutations_logreg)
    df_importances_logreg.to_parquet(args.out_file_df_importances_logreg)
    df_markers_rank_significance_logreg.to_parquet(args.out_file_df_markers_rank_significance_logreg)

    # diffexpr
    df_diffexpr.to_parquet(args.out_file_df_diffexpr)
    df_markers_rank_significance_diffexpr.to_parquet(args.out_file_df_markers_rank_significance_diffexpr)

    # liana
    # readwrite.write_anndata_folder(lrdata, args.out_dir_liana_lrdata)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

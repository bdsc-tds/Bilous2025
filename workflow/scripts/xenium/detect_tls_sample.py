import dask

dask.config.set({"dataframe.query-planning": False})

import numpy as np
import pandas as pd
import argparse
import os
import sys

sys.path.append("workflow/scripts/")
import tls
import readwrite


# Set up argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Run RESOLVI on a Xenium sample.")
    parser.add_argument("--path", type=str, help="Path to the xenium sample file.")
    parser.add_argument("--out_file", type=str, help="output file")
    parser.add_argument("--radius", type=float, help="")
    parser.add_argument("--b_cell_label", type=str, help="")
    parser.add_argument("--t_cell_label", type=str, help="")
    parser.add_argument("--min_perc_b", type=float, help="")
    parser.add_argument("--min_perc_t", type=float, help="")
    parser.add_argument("--fold_change_threshold", type=float, help="")
    parser.add_argument("--min_b", type=int, help="")
    parser.add_argument("--min_t", type=int, help="")
    parser.add_argument("--min_tls_size", type=int, help="")
    parser.add_argument("--cell_type_labels", type=str, help="cell_type_labels")
    parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
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

    adata.X.data = adata.X.data.astype(np.float32)
    adata.obs_names = adata.obs_names.astype(str)

    if "proseg" in args.path and "raw_results" in args.path:
        # need to round proseg expected counts for resolVI to run
        adata.X.data = adata.X.data.round()
        adata.obs_names = "proseg-" + adata.obs_names

    label_key = "cell_type"
    adata.obs[label_key] = pd.read_parquet(args.cell_type_labels).set_index("cell_id").iloc[:, 0].astype("str")
    adata = adata[adata.obs[label_key].notna()]

    if "Level2.1" in args.cell_type_labels:
        # for custom Level2.1, simplify subtypes
        adata.obs.loc[adata.obs[label_key].str.contains("malignant"), label_key] = "malignant cell"
        adata.obs.loc[adata.obs[label_key].str.contains("T cell"), label_key] = "T cell"

    tls.detect_tls(
        adata,
        radius=args.radius,
        b_cell_label=args.b_cell_label,
        t_cell_label=args.t_cell_label,
        min_perc_b=args.min_perc_b,
        min_perc_t=args.min_perc_t,
        fold_change_threshold=args.fold_change_threshold,
        min_b=args.min_b,
        min_t=args.min_t,
        min_tls_size=args.min_tls_size,
        cell_type_obs_key=label_key,
        spatial_obsm_key="spatial",
        output_obs_key="tls_cluster",
    )

    columns = ["neigh_count_b", "neigh_count_t", "neigh_prop_b", "neigh_prop_t", "is_candidate", "tls_cluster"]
    df = adata.obs[columns]
    df.to_parquet(args.out_file)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

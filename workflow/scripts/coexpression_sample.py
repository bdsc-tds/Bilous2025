"""
Compute coexpression with given method and target counts per segmentation method per sample.
"""

import dask

dask.config.set({"dataframe.query-planning": False})

import os
import argparse

from pathlib import Path
import scanpy as sc
import sys

sys.path.append("workflow/scripts/")
import coexpression as ce
import readwrite


def parse_args():
    sys_args_parser = argparse.ArgumentParser(description="Compute coexpression.")

    sys_args_parser.add_argument(
        "-i",
        required=True,
        type=Path,
        help="path to the 10X XeniumRanger output folder",
    )
    sys_args_parser.add_argument(
        "-l",
        type=str,
        default=None,
        help="path to the log file",
    )
    sys_args_parser.add_argument(
        "--outcoexp",
        required=True,
        type=str,
        help="path to the output file for coexpression",
    )
    sys_args_parser.add_argument(
        "--outposrate",
        required=True,
        type=str,
        help="path to the output file for positive rate",
    )
    sys_args_parser.add_argument(
        "-m",
        required=True,
        type=str,
        help="method to compute coexpression",
    )
    sys_args_parser.add_argument(
        "-c",
        required=True,
        type=int,
        help="target count for coexpression",
    )
    sys_args_parser.add_argument(
        "-g",
        required=False,
        type=str,
        help="optional genes subset to use: 'condition' (common across condition panels) or 'conditions' (common across all conditions' panels). ",
    )

    ret = sys_args_parser.parse_args()
    if not os.path.isdir(ret.i) and ret.i.suffix != ".h5":
        raise RuntimeError(f"Error! 10X XeniumRanger output folder does not exist: {ret.i}")
    if ret.i.suffix == ".h5" and not ret.i.exists():
        raise RuntimeError(f"Error! 10X h5 file does not exist: {ret.i}")
    os.makedirs(os.path.dirname(ret.outcoexp), exist_ok=True)
    os.makedirs(os.path.dirname(ret.outposrate), exist_ok=True)

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
    if args.i.suffix == ".h5":
        adata = sc.read_10x_h5(args.i)

    else:
        adata = readwrite.read_xenium_sample(
            args.i,
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

    # prevent scanpy bug when subsetting with int index
    adata.obs_names = adata.obs_names.astype(str)

    # subset to common genes across panels for a condition or across all conditions
    if args.g in ["condition", "conditions"]:
        if "raw_results" in args.i.as_posix():
            # proseg expected counts directory
            xenium_processed_data_dir = args.i.parents[5]
            sample_condition = args.i.parents[4].stem
        else:
            xenium_processed_data_dir = args.i.parents[6]
            sample_condition = args.i.parents[5].stem

        ref_segmentation = list(xenium_processed_data_dir.iterdir())[0]  # get any segmentation

        panels_genes = {}
        for condition in ref_segmentation.iterdir():
            # only use panels within the same condition if args.g == "condition"
            if args.g == "condition" and condition.stem != sample_condition:
                continue
            for panel in condition.iterdir():
                # get any sample
                donor = list(panel.iterdir())[0]
                sample = list(donor.iterdir())[0]

                name = "/".join((condition.stem, panel.stem, donor.stem, sample.stem))
                p = ref_segmentation / f"{name}/normalised_results/outs/gene_panel.json"
                panels_genes[condition.stem, panel.stem] = readwrite.get_gene_panel_info(p).query("id.notnull()")[
                    "name"
                ]
        common_genes = list(set.intersection(*map(set, panels_genes.values())))

    elif isinstance(args.g, str):
        # read gene panel names
        gene_panel = readwrite.get_gene_panel_info(args.g)["name"]
        common_genes = [gene for gene in gene_panel if gene in adata.var_names]
    else:
        common_genes = None

    # subset adata
    if common_genes is not None:
        print(adata.obs_names)
        print(type(adata.obs_names))

        adata = adata[:, common_genes].copy()
        if adata.shape[1] == 0:
            raise ValueError("No genes found from adata.var_names")

    # compute coexpression
    CC, X_downsampled, pos, pos_rate, mask = ce.coexpression(
        adata,
        target_count=args.c,
        method=args.m,
    )

    # save as parquet
    CC.to_parquet(args.outcoexp)
    pos_rate.to_frame().to_parquet(args.outposrate)

    if args.l is not None:
        _log.close()
        sys.stdout = old_stdout
        sys.stderr = old_stderr

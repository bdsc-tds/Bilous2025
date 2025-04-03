# %%
import dask

dask.config.set({"dataframe.query-planning": False})

import argparse
import numpy as np
import pandas as pd
import sys
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append("workflow/scripts/")
import _utils
import readwrite

sns.set_style("ticks")

# Set up argument parser
parser = argparse.ArgumentParser(description="")
parser.add_argument("--condition", type=str, help="condition name.")
parser.add_argument("--panel", type=str, help="panel name.")
parser.add_argument("--correction_methods", type=str, nargs="*", default=[], help="correction methods list.")
parser.add_argument("--results_dir", type=Path, help="Path to the results dir.")
parser.add_argument("--std_seurat_analysis_dir", type=Path, help="Path to the std seurat analysis dir.")
parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell type annotation dir.")
parser.add_argument("--out_dir", type=Path, help="Path to the output dir.")
parser.add_argument("--normalisation", type=str, help="Normalisation method")
parser.add_argument("--layer", type=str, help="Name of saved layer of the seurat object, data or scale_data")
parser.add_argument("--reference", type=str, help="annotation reference")
parser.add_argument("--method", type=str, help="annotation method")
parser.add_argument("--level", type=str, help="annotation level")
parser.add_argument("--n_comps", type=int, help="Number of components.")
parser.add_argument("--max_n_cells", type=int, help="Max number of cells to use.")
parser.add_argument("--top_n", type=int, help="contamination diffexpr script parameter")
parser.add_argument("--mixture_k", type=int, help="ResolVI parameter")
parser.add_argument("--num_samples", type=int, help="ResolVI parameter")
parser.add_argument("--use_precomputed", action="store_true", help="Use precomputed data.")
parser.add_argument("--count_correction_palette", type=Path, help="Path to the count correction palette file.")
parser.add_argument("--dpi", type=int, help="Figure DPI.")
parser.add_argument("--extension", type=str, help="conservation metric to plot")
args = parser.parse_args()

# Access the arguments
condition = args.condition
panel = args.panel
correction_methods = args.correction_methods
results_dir = args.results_dir
cell_type_annotation_dir = args.cell_type_annotation_dir
std_seurat_analysis_dir = args.std_seurat_analysis_dir
out_dir = args.out_dir
normalisation = args.normalisation
layer = args.layer
reference = args.reference
method = args.method
level = args.level
n_comps = args.n_comps
max_n_cells = args.max_n_cells
top_n = args.top_n
mixture_k = args.mixture_k
num_samples = args.num_samples
count_correction_palette = args.count_correction_palette
use_precomputed = args.use_precomputed
dpi = args.dpi
extension = args.extension
args = parser.parse_args()


# Params
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
order = ["breast", "chuvio", "lung", "5k"]
hue_segmentation = "segmentation"
hue_segmentation_order = [
    "10x_mm_0um",
    "10x_mm_5um",
    "10x_mm_15um",
    "10x_0um",
    "10x_5um",
    "10x_15um",
    "baysor",
    "proseg_expected",
    "proseg_mode",
    "segger",
]
hue_correction = "correction_method"
hue_correction_order = [
    "raw",
    "ResolVI",
    "ResolVI supervised",
    "ovrlpy 0.5",
    "ovrlpy 0.7",
    "SPLIT",
]

rank_metrics = ["logfoldchanges", "-log10pvals_x_logfoldchanges", "-log10pvals_x_sign_logfoldchanges", "mean_zscore"]
plot_metrics = ["n_cells", "mean_n_counts", "mean_n_genes", "median_n_genes"]

# %% [markdown]
# # Load results diffexpr

# %%
dfs = readwrite.read_contamination_metrics_results(
    results_dir,
    correction_methods,
    std_seurat_analysis_dir,
    reference,
    method,
    level,
    mixture_k,
    num_samples,
    normalisation,
    layer,
    ref_condition=condition,
    ref_panel=panel,
    evaluation="diffexpr",
)


# %% [markdown]
# # Plot decontamination results diffexpr

# %%

df_count_correction_palette = pd.read_csv(count_correction_palette, index_col=0).iloc[:, 0]


for plot_metric in plot_metrics:
    out_file = out_dir / f"{panel}_{plot_metric}.png"

    df = _utils.get_df_summary_stats_plot(dfs, plot_metric=plot_metric)
    df = df.query("panel == @panel")

    # rename proseg
    df.loc[df["segmentation"] == "proseg_expected", "segmentation"] = "proseg"
    df_count_correction_palette = df_count_correction_palette.rename(index={"proseg_expected": "proseg"})

    # plotting params, palette
    unique_labels = [c for c in hue_correction_order if c in np.unique(df[hue_correction].dropna())]
    unique_labels = unique_labels + [c for c in np.unique(df[hue_correction].dropna()) if c not in unique_labels]
    palette = {u: df_count_correction_palette[u] for u in unique_labels}
    legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]

    ###  boxplot
    f = plt.figure(figsize=(5, 3))
    ax = plt.subplot()
    g = sns.boxplot(
        df,
        x="segmentation",
        y=plot_metric,
        hue=hue_correction,
        hue_order=unique_labels,
        legend=False,
        palette=palette,
        ax=ax,
        order=[s for s in hue_segmentation_order if s in df["segmentation"].unique()],
        log_scale=True if plot_metric == "n_cells" else False,
        showfliers=True,
    )

    if plot_metric == "n_cells":
        ax.set_ylim(None, 1e6)
    sns.despine(offset=10, trim=True)
    ax.yaxis.grid(True)
    ax.yaxis.set_tick_params(labelsize=12)  # If you also want to change the y-axis numbers
    if plot_metric == '"-log10pvalue"':
        ax.set.ylabel(r"$-\log_{10} \text{ p-value}$", fontsize=14)
    else:
        ax.set_ylabel(plot_metric, fontsize=14)
    plt.setp(ax.get_xticklabels(), rotation=45, fontsize=12)

    # title = f"{plot_metric} for {panel}"
    # plt.suptitle(title)
    # f.legend(
    #     handles=legend_handles,
    #     loc="center left",
    #     bbox_to_anchor=(1, 0.5),
    #     title=hue_correction,
    #     frameon=False,
    # )
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight")

from pathlib import Path
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import natsort

sys.path.append("workflow/scripts/")
import readwrite
import coexpression


def float_or_str(value):
    if value == "auto":
        return value
    try:
        return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid value: {value}. Must be a float or 'auto'.")


def format_ticks(x):
    if x < 1:
        if "5" in str(x.round(6)) or "1" in str(x.round(6)):
            return f".{str(x.round(6))[2:]}"
        else:
            return ""
    else:
        if "5" in str(x.round(6)) or "1" in str(x.round(6)):
            return int(x)
        else:
            return ""


# Set up argument parser
parser = argparse.ArgumentParser(description="Plot violins of Xenium coexpression QC for a given panel.")
parser.add_argument("--coexpression_dir", type=Path, help="Path to coexpression results to plot")
parser.add_argument("--plot_condition", type=str, help="condition results to plot")
parser.add_argument("--plot_panel", type=str, help=" panel results to plot")
parser.add_argument("--out_file_plot_sample", type=Path, help="Path to the output plot file")
parser.add_argument("--out_file_plot_panel", type=Path, help="Path to the output plot file")
parser.add_argument("--out_file_gene_pairs", type=Path, help="Path to the output gene pairs file")
parser.add_argument("--method", type=str, help="method annotation parameter")
parser.add_argument("--target_count", type=int, help="target count parameter")
parser.add_argument("--min_positivity_rate", type=float, help="min_positivity_rate")
parser.add_argument("--min_cond_coex", type=float_or_str, help="min_cond_coex")
parser.add_argument("--min_cond_coex_mode", type=str, help="min_cond_coex_mode")
parser.add_argument("--cc_cutoff", type=float, help="coexpression cutoff")
parser.add_argument("--ref_segmentation", type=str, help="reference segmentation")
parser.add_argument("--ref_oversegmentation", type=str, help="reference oversegmentation")
parser.add_argument("--segmentation_palette", type=Path, help="segmentation palette")
parser.add_argument("--dpi", type=int, help="figures dpi")
parser.add_argument("--showfliers", action="store_true", help="showfliers or not in boxplots")
parser.add_argument("--log_scale", action="store_true", help="log_scale for boxplots")

args = parser.parse_args()

# Access the arguments
coexpression_dir = args.coexpression_dir
plot_condition = args.plot_condition
plot_panel = args.plot_panel
out_file_plot_sample = args.out_file_plot_sample
out_file_plot_panel = args.out_file_plot_panel
out_file_gene_pairs = args.out_file_gene_pairs
method = args.method
target_count = args.target_count
min_positivity_rate = args.min_positivity_rate
min_cond_coex = args.min_cond_coex
min_cond_coex_mode = args.min_cond_coex_mode
cc_cutoff = args.cc_cutoff
ref_segmentation = args.ref_segmentation
ref_oversegmentation = args.ref_oversegmentation
palette = pd.read_csv(args.segmentation_palette, index_col=0).iloc[:, 0]
dpi = args.dpi
showfliers = args.showfliers
log_scale = args.log_scale

print(f"{showfliers=}")

# vars
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
hue = "segmentation"
hue_order = ["10x_mm_0um", "10x_mm_5um", "10x_mm_15um", "10x_0um", "10x_5um", "10x_15um", "baysor", "proseg", "segger"]


# Read coexpression results
cc_paths = []
for segmentation in (segmentations := coexpression_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()):
        if condition.stem != plot_condition:
            continue
        for panel in (panels := condition.iterdir()):
            if panel.stem != plot_panel:
                continue
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    k = (
                        segmentation.stem,
                        condition.stem,
                        panel.stem,
                        donor.stem,
                        sample.stem,
                    )

                    cc_paths.append((k, method, target_count))

CC, pos_rate = readwrite.read_coexpression_files(cc_paths, coexpression_dir)


# automatically select min_cond_coex
if method in ["conditional", "jaccard"] and min_cond_coex == "auto":
    min_cond_coex = float(np.min([np.nanmin(CC[k][method, target_count].replace(0.0, np.nan)) for k in CC.keys()]))


# compute ratio with ref_segmentation
CCdiff = {}
spurious_gene_pairs = {}
for k in CC.keys():
    if k[0] == ref_segmentation:
        continue

    CCdiff[k] = {}
    spurious_gene_pairs[k] = {}

    k_ref = (ref_segmentation, *k[1:])

    (CCdiff[k][method, target_count], spurious_gene_pairs[k][method, target_count]) = (
        coexpression.compare_segmentations(
            CC_ref_seg=CC[k_ref][method, target_count],
            CC_other_seg=CC[k][method, target_count],
            pos_rate_ref_seg=pos_rate[k_ref][method, target_count],
            pos_rate_other_seg=pos_rate[k][method, target_count],
            min_positivity_rate=min_positivity_rate,
            min_cond_coex=min_cond_coex,
            min_cond_coex_mode=min_cond_coex_mode,
            cc_cutoff=cc_cutoff,
            method=method,
        )
    )


# extract spurious gene pairs based on ref_oversegmentation/ref_segmentation ratio
# and stack into one df
keys = pd.DataFrame(CCdiff.keys(), columns=xenium_levels)

data = []
for _, k in keys.iterrows():
    k_ref_over = (ref_oversegmentation, *k[1:])

    mat = CCdiff[*k][method, target_count]

    i = spurious_gene_pairs[k_ref_over][method, target_count][:, 0]
    j = spurious_gene_pairs[k_ref_over][method, target_count][:, 1]
    flat_values = mat.values[mat.index.get_indexer(i), mat.columns.get_indexer(j)]
    data.extend(np.hstack((np.tile(k, (len(i), 1)), i[:, None], j[:, None], flat_values[:, None])))


# Convert to DataFrame for plotting
df = pd.DataFrame(data, columns=xenium_levels + ["genei", "genej", "relative coexpression"])
df["relative coexpression"] = df["relative coexpression"].astype(float)


# plotting params, palette
unique_labels = [ct for ct in hue_order if ct in np.unique(df[hue].dropna())]
unique_labels = unique_labels + [ct for ct in np.unique(df[hue].dropna()) if ct not in unique_labels]
palette = {u: palette[u] for u in unique_labels}


# Create boxplot
for y in ["sample", "panel"]:
    if y == "sample":
        figsize = (6, df["sample"].nunique())
        out_file_ = out_file_plot_sample
    else:
        figsize = (6, 6)
        out_file_ = out_file_plot_panel

    print(f"{showfliers=}")

    plt.figure(figsize=figsize)
    ax = plt.subplot()
    g = sns.boxplot(
        data=df,
        x="relative coexpression",
        y=y,
        hue=hue,
        hue_order=unique_labels,
        palette=palette,
        # cut=0,
        # width=0.8,
        # inner="quart",
        log_scale=log_scale,
        showfliers=showfliers,
        flierprops={
            "marker": "o",
            "color": "black",
            "markersize": 1,
            "markerfacecolor": "w",
        },
        ax=ax,
    )

    # sns.despine(offset=10, trim=True)

    if log_scale:
        ax.set_xticklabels([format_ticks(x) for x in ax.get_xticks(minor=True)], minor=True)
        ax.set_xticklabels([format_ticks(x) for x in ax.get_xticks()])
        ax.tick_params(axis="both", which="major", labelsize=14)
        ax.tick_params(axis="both", which="minor", labelsize=10)

    ax.xaxis.grid(True)
    plt.title(f"{plot_condition=} {plot_panel=} {method=} {target_count=}")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), title=hue, frameon=False)
    plt.savefig(out_file_, dpi=dpi, bbox_inches="tight")
    plt.close()

df.to_parquet(out_file_gene_pairs)

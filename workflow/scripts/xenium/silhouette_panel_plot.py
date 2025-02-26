import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Set up argument parser
parser = argparse.ArgumentParser(description="Plot silhouette of Xenium donors for a given panel.")
parser.add_argument("--silhouette_dir", type=Path, help="Path to silhouette_dir to plot")
parser.add_argument("--plot_panel", type=str, help="name of panel to plot")
parser.add_argument("--plot_condition", type=str, help="name of condition to plot")
parser.add_argument("--out_file", type=Path, help="Path to the output file")
parser.add_argument("--segmentation_palette", type=Path, help="segmentation palette")
parser.add_argument("--normalisation", type=str, help="normalisation method")
parser.add_argument("--layer", type=str, help="layer method")
parser.add_argument("--reference", type=Path, help="reference annotation parameter")
parser.add_argument("--method", type=Path, help="method annotation parameter")
parser.add_argument("--level", type=Path, help="level annotation parameter")
parser.add_argument("--dpi", type=int, help="dpi of saved plot")
parser.add_argument("--metric", type=int, help="distance metric used to compute silhouettes")

args = parser.parse_args()

# Access the arguments
silhouette_dir = args.silhouette_dir
plot_condition = args.plot_condition
plot_panel = args.plot_panel
out_file = args.out_file
segmentation_palette = args.segmentation_palette
normalisation = args.normalisation
layer = args.layer
reference = args.reference
method = args.method
level = args.level
dpi = args.dpi
metric = args.metric

# vars
palette = pd.read_csv(segmentation_palette, index_col=0)["cols_segmentation"]
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample", "index"]
hue = "segmentation"
hue_order = [
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

# read cell type annotation
annot = {}
for segmentation in (segmentations := silhouette_dir.iterdir()):
    print(segmentation.stem)
    for condition in (conditions := segmentation.iterdir()):
        if condition.stem != plot_condition:
            continue
        for panel in (panels := condition.iterdir()):
            if panel.stem != plot_panel:
                continue
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    k = (segmentation.stem, condition.stem, panel.stem, donor.stem, sample.stem)

                    annot[k] = {}
                    annot_file = sample / f"{normalisation}/silhouette_{layer}_{metric}.parquet"

                    # if annot_file.exists():
                    annot[k] = pd.read_parquet(annot_file)


# merge annotations
df_annot = pd.concat(annot)
df_annot = df_annot.reset_index()
df_annot.columns = xenium_levels + df_annot.columns[len(xenium_levels) :].tolist()

df = df_annot.query(f"reference == '{reference}' and method == '{method}' and level == '{level}'")

# plotting params, palette
title = f"condition: {plot_condition}, Panel: {plot_panel}\n Reference: {reference}, Method: {method}, Level: {level}"
unique_labels = [c for c in hue_order if c in np.unique(df[hue].dropna())]
unique_labels = unique_labels + [c for c in np.unique(df[hue].dropna()) if c not in unique_labels]
palette = {u: palette[u] for u in unique_labels}
legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]


# Create joint boxplot
sns.set(style="ticks")
f = plt.figure(figsize=(6, df["sample"].nunique() // 2))
g = sns.stripplot(data=df, x="silhouette", y="sample", hue=hue, hue_order=unique_labels, legend=False, palette=palette)

sns.despine(offset=10, trim=True)
plt.gca().xaxis.grid(True)
plt.axvline(0, c="k", linestyle="-", alpha=0.4)

plt.title(title)
f.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    title=hue,
    frameon=False,
)
plt.tight_layout(rect=[0, 0, 0.85, 0.95])
plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
plt.close()

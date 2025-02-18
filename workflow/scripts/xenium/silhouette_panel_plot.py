import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Set up argument parser
parser = argparse.ArgumentParser(description="Plot silhouette of Xenium donors for a given panel.")
parser.add_argument("--panel", type=Path, help="Path to panel to plot")
parser.add_argument("--out_file", type=Path, help="Path to the output file")
parser.add_argument("--segmentation_palette", type=Path, help="segmentation palette")
parser.add_argument("--normalisation", type=str, help="normalisation method")
parser.add_argument("--layer", type=str, help="layer method")
parser.add_argument("--reference", type=Path, help="reference annotation parameter")
parser.add_argument("--method", type=Path, help="method annotation parameter")
parser.add_argument("--level", type=Path, help="level annotation parameter")
args = parser.parse_args()

# Access the arguments
panel = args.panel
out_file = args.out_file
segmentation_palette = args.segmentation_palette
normalisation = args.normalisation
layer = args.layer
reference = args.reference
method = args.method
level = args.level

# vars
silhouette_dir = panel.parents[2]
plot_panel = panel.stem
plot_condition = panel.parents[0].stem
palette = pd.read_csv(segmentation_palette, index_col=0)
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample", "index"]
hue = "segmentation"

# read cell type annotation
annot = {}
for segmentation in (segmentations := silhouette_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()):
        if condition.stem != plot_condition:
            continue
        for panel in (panels := condition.iterdir()):
            if panel.stem != plot_panel:
                continue
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    k = (segmentation.stem, condition.stem, panel.stem, donor.stem, sample.stem, normalisation)

                    annot[k] = {}
                    annot_file = sample / f"silhouette_{layer}.parquet"
                    if annot_file.exists():
                        annot[k] = pd.read_parquet(annot_file)

# merge annotations
df_annot = pd.concat(annot)
df_annot = df_annot.reset_index()
df_annot.columns = xenium_levels + df_annot.columns[len(xenium_levels) :].tolist()

df = df_annot.query(  # condition == '{plot_condition}' and panel == '{plot_panel}' and
    f"reference == '{reference}' and method == '{method}' and level == '{level}'"
)

# plotting params, palette
title = f"condition: {plot_condition}, Panel: {plot_panel}\n Reference: {reference}, Method: {method}, Level: {level}"
unique_labels = np.unique(df[hue].dropna())
palette = {u: palette[u] for u in unique_labels}
legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]


# Create joint boxplot
sns.set(style="ticks")
f = plt.figure(figsize=(6, df["sample"].nunique() // 2))
g = sns.stripplot(data=df, x="silhouette", y="sample", hue=hue, legend=False, palette=palette)

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
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.close()

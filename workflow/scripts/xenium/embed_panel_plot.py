import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Set up argument parser
parser = argparse.ArgumentParser(description="Plot panel of Xenium samples.")
parser.add_argument("--panel", type=Path, help="Path to the panel file.")
parser.add_argument("--embed_file", type=str, help="Path to the embedding file.")
parser.add_argument("--reference", type=str, help="annotation reference")
parser.add_argument("--method", type=str, help="annotation method")
parser.add_argument("--level", type=str, help="annotation level")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
args = parser.parse_args()

# Access the arguments
panel = args.panel
embed_file = args.embed_file
reference = args.reference
method = args.method
level = args.level
out_file = args.out_file

# vars
xenium_levels = ["segmentation", "cohort", "panel", "sample", "replicate", "cell_id"]
segmentation = panel.parents[1]
cohort = panel.parents[0]

# load umap
obs = pd.read_parquet(embed_file)
obs["cell_id"] = obs.index


if level == "replicate":
    # plot replicate as color, no need to load annotations
    df = obs
    params = level
    title = (
        f"Segmentation: {segmentation.stem}, Cohort: {cohort.stem}, Panel: {panel.stem}"
    )

else:
    # read cell type annotation
    annot = {}
    for sample in (samples := panel.iterdir()):
        for replicate in (replicates := sample.iterdir()):
            k = (
                segmentation.stem,
                cohort.stem,
                panel.stem,
                sample.stem,
                replicate.stem,
            )
            name = "/".join(k)

            annot[k] = {}
            annot_file = (
                replicate
                / f"cell_type_annotation/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet"
            )
            if annot_file.exists():
                annot[k][reference, method, level] = pd.read_parquet(annot_file).iloc[
                    :, 0
                ]

            annot_file = (
                replicate
                / f"cell_type_annotation/reference_based/{reference}/{method}/{level}/single_cell/labels.csv"
            )
            if annot_file.exists():
                annot[k][reference, method, level] = pd.read_csv(
                    annot_file, index_col=0
                ).iloc[:, 0]

    # merge annotations
    df_annot = {}
    for k in annot:
        if len(annot[k]):
            df_ = pd.DataFrame(annot[k])
            df_.columns = [col for col in df_.columns]
            df_annot[k] = df_
    df_annot = pd.concat(df_annot)
    df_annot.index.names = xenium_levels
    df_annot = df_annot.reset_index()

    # merge umap and cell type annotations
    df = pd.merge(obs, df_annot, on=xenium_levels, how="inner")

    params = (reference, method, level)
    title = f"Segmentation: {segmentation.stem}, Cohort: {cohort.stem}, Panel: {panel.stem}\n Method: {method}, Reference: {reference}, Level: {level}"


# plotting params, palette
unique_labels = np.unique(df[params].dropna())
palette = dict(zip(unique_labels, sc.pl.palettes.default_28))
legend_handles = [
    mpatches.Patch(color=color, label=label) for label, color in palette.items()
]

print(
    f"Segmentation: {segmentation.stem}, Cohort: {cohort.stem}, Panel: {panel.stem}, Method: {method}, Reference: {reference}, Level: {level}"
)


# plot
f = plt.figure(figsize=(12, 10))
ax = plt.subplot()

sns.scatterplot(
    data=df,
    x="UMAP1",
    y="UMAP2",
    s=0.1,
    alpha=0.5,
    hue=params,
    ax=ax,
    palette=palette,
    legend=False,
)
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
sns.despine()

plt.title(title)
f.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    title=params if isinstance(params, str) else ", ".join(params),
    frameon=False,
)
plt.tight_layout(rect=[0, 0, 0.85, 0.95])
plt.savefig(out_file, dpi=300, bbox_inches="tight")
plt.close()

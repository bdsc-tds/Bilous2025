import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Set up argument parser
parser = argparse.ArgumentParser(description="Plot condition of Xenium donors.")
parser.add_argument("--condition", type=Path, help="Path to the condition file.")
parser.add_argument("--embed_file", type=str, help="Path to the embedding file.")
parser.add_argument("--cell_type_annotation_dir", type=Path, help="Path to the cell_type_annotation_dir.")
parser.add_argument("--normalisation", type=str, help="normalisation method")
parser.add_argument("--reference", type=str, help="annotation reference")
parser.add_argument("--method", type=str, help="annotation method")
parser.add_argument("--color", type=str, help="annotation color")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--cell_type_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--panel_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--sample_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--s", type=float, help="scatter point size")
parser.add_argument("--alpha", type=float, help="scatter alpha (transparency)")
args = parser.parse_args()

# Access the arguments
condition = args.condition
embed_file = args.embed_file
cell_type_annotation_dir = args.cell_type_annotation_dir
normalisation = args.normalisation
reference = args.reference
method = args.method
color = args.color
out_file = args.out_file
cell_type_palette = args.cell_type_palette
panel_palette = args.panel_palette
sample_palette = args.sample_palette
s = args.s
alpha = args.alpha

if color == "sample":
    palette = pd.read_csv(sample_palette, index_col=0).iloc[:, 0]
elif color == "panel":
    palette = pd.read_csv(panel_palette, index_col=0).iloc[:, 0]
else:
    palette = pd.read_csv(cell_type_palette)[[color, f"cols_{color}"]].drop_duplicates().set_index(color).squeeze()


# vars
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample", "cell_id"]
segmentation = condition.parents[0]

# load umap
obs = pd.read_parquet(embed_file)
obs["cell_id"] = obs.index


if color in ["panel", "sample"]:
    # plot sample as color, no need to load annotations
    df = obs
    params = color
    title = f"Segmentation: {segmentation.stem}, condition: {condition.stem}"

else:
    # read cell type annotation
    annot = {}

    for panel in (panels := condition.iterdir()):
        for donor in (donors := panel.iterdir()):
            for sample in (samples := donor.iterdir()):
                k = (
                    segmentation.stem,
                    condition.stem,
                    panel.stem,
                    donor.stem,
                    sample.stem,
                )
                name = "/".join(k)

                annot[k] = {}
                annot_file = (
                    cell_type_annotation_dir
                    / name
                    / f"{normalisation}/reference_based/{reference}/{method}/{color}/single_cell/labels.parquet"
                )
                # if annot_file.exists():
                annot[k][reference, method, color] = pd.read_parquet(annot_file).set_index("cell_id").iloc[:, 0]

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

    params = (reference, method, color)
    title = f"Segmentation: {segmentation.stem}, condition: {condition.stem}\n Method: {method}, Reference: {reference}"


# plotting params, palette
unique_labels = np.unique(df[params].dropna())
palette = {u: palette[u] for u in unique_labels}
legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]

print(f"Segmentation: {segmentation.stem}, condition: {condition.stem}, Method: {method}, Reference: {reference}")


# plot
f = plt.figure(figsize=(12, 10))
ax = plt.subplot()

sns.scatterplot(
    data=df,
    x="UMAP1",
    y="UMAP2",
    s=s,
    alpha=alpha,
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

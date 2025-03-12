import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# Set up argument parser
parser = argparse.ArgumentParser(description="Plot panel of Xenium donors.")
parser.add_argument("--embed_file", type=str, help="Path to the embedding file.")
parser.add_argument("--reference", type=Path, help="annotation reference")
parser.add_argument("--color", type=str, help="annotation color")
parser.add_argument("--out_file", type=str, help="Path to the output file.")
parser.add_argument("--cell_type_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--panel_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--sample_palette", type=Path, help="Path to palette csv file")
parser.add_argument("--s", type=float, help="scatter point size")
parser.add_argument("--alpha", type=float, help="scatter alpha (transparency)")
parser.add_argument("--dpi", type=int, help="dpi of saved plot")
args = parser.parse_args()

# Access the arguments
embed_file = args.embed_file
reference = args.reference
color = args.color
out_file = args.out_file
cell_type_palette = args.cell_type_palette
panel_palette = args.panel_palette
sample_palette = args.sample_palette
s = args.s
alpha = args.alpha
dpi = args.dpi

if color == "sample":
    palette = pd.read_csv(sample_palette, index_col=0).iloc[:, 0]
elif color == "panel":
    palette = pd.read_csv(panel_palette, index_col=0).iloc[:, 0]
else:
    palette = pd.read_csv(cell_type_palette)[[color, f"cols_{color}"]].drop_duplicates().set_index(color).squeeze()


# load umap
df = pd.read_parquet(embed_file)
df["cell_id"] = df.index
df[color] = pd.read_parquet(reference / "metadata.parquet")[color].values

# plotting color, palette
unique_labels = np.unique(df[color].dropna())
palette = {u: palette[u] for u in unique_labels}
legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]

# plot
f = plt.figure(figsize=(12, 10))
ax = plt.subplot()

sns.scatterplot(
    data=df,
    x="UMAP1",
    y="UMAP2",
    s=s,
    alpha=alpha,
    hue=color,
    ax=ax,
    palette=palette,
    legend=False,
)
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
sns.despine()

f.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    title=color if isinstance(color, str) else ", ".join(color),
    frameon=False,
)
plt.tight_layout(rect=[0, 0, 0.85, 0.95])
plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
plt.close()

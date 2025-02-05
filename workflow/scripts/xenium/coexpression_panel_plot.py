from pathlib import Path
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

sys.path.append("workflow/scripts/")
import readwrite
import coexpression

# Set up argument parser
parser = argparse.ArgumentParser(
    description="Plot violins of Xenium coexpression QC for a given panel."
)
parser.add_argument("--panel", type=Path, help="Path to panel results to plot")
parser.add_argument(
    "--out_file_plot_sample", type=Path, help="Path to the output plot file"
)
parser.add_argument(
    "--out_file_plot_panel", type=Path, help="Path to the output plot file"
)
parser.add_argument(
    "--out_file_gene_pairs", type=Path, help="Path to the output gene pairs file"
)
parser.add_argument("--method", type=str, help="method annotation parameter")
parser.add_argument("--target_count", type=int, help="target count parameter")
parser.add_argument("--min_positivity_rate", type=float, help="min_positivity_rate")
# parser.add_argument("--min_cond_coex", type=float, help="min_cond_coex")
parser.add_argument("--cc_cutoff", type=float, help="coexpression cutoff")
parser.add_argument("--log2", type=bool, help="log2 scale or not")
parser.add_argument("--ref_segmentation", type=str, help="reference segmentation")
parser.add_argument(
    "--ref_oversegmentation", type=str, help="reference oversegmentation"
)
parser.add_argument("--segmentation_palette", type=Path, help="segmentation palette")

args = parser.parse_args()

# Access the arguments
panel = args.panel
out_file_plot_sample = args.out_file_plot_sample
out_file_plot_panel = args.out_file_plot_panel
out_file_gene_pairs = args.out_file_gene_pairs
method = args.method
target_count = args.target_count
min_positivity_rate = args.min_positivity_rate
# min_cond_coex = args.min_cond_coex
cc_cutoff = args.cc_cutoff
log2 = args.log2
ref_segmentation = args.ref_segmentation
ref_oversegmentation = args.ref_oversegmentation

palette = pd.read_csv(args.segmentation_palette, index_col=0).iloc[:, 0]

# vars
coexpression_dir = panel.parents[2]
plot_panel = panel.stem
plot_condition = panel.parents[0].stem
xenium_levels = ["segmentation", "condition", "panel", "donor", "sample"]
hue = "segmentation"

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

if method in ["conditional", "jaccard"]:
    min_cond_coex = float(
        np.min(
            [
                np.nanmin(CC[k][method, target_count].replace(0.0, np.nan))
                for k in CC.keys()
            ]
        )
    )
else:
    min_cond_coex = 0.0

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
            cc_cutoff=cc_cutoff,
            method=method,
            log2=False,
        )
    )


# extract spurious gene pairs based on ref_oversegmentation/ref_segmentation ratio
# and stack into one df
keys = pd.DataFrame(CCdiff.keys(), columns=xenium_levels)

data = []
for _, k in keys.iterrows():
    k_ref_over = (ref_oversegmentation, *k[1:])

    mat = CCdiff[*k][method, target_count]
    if log2:
        with warnings.catch_warnings(action="ignore"):
            mat = np.log2(mat)

    i = spurious_gene_pairs[k_ref_over][method, target_count][:, 0]
    j = spurious_gene_pairs[k_ref_over][method, target_count][:, 1]
    flat_values = mat.values[mat.index.get_indexer(i), mat.columns.get_indexer(j)]
    data.extend(
        np.hstack(
            (np.tile(k, (len(i), 1)), i[:, None], j[:, None], flat_values[:, None])
        )
    )


# Convert to DataFrame for plotting
df = pd.DataFrame(
    data, columns=xenium_levels + ["genei", "genej", "log2 relative coexpression"]
)
df["log2 relative coexpression"] = df["log2 relative coexpression"].astype(float)

# set inf values to max
# max = df["log2 relative coexpression"][df["log2 relative coexpression"]!=np.inf].max()
# df["log2 relative coexpression"] = df["log2 relative coexpression"].replace(np.inf, max)


# plotting params, palette
unique_labels = np.unique(df[hue].dropna())
palette = {u: palette[u] for u in unique_labels}

# Create sample violinplot
plt.figure(figsize=(6, df["sample"].nunique()))
g = sns.violinplot(
    data=df,
    x="log2 relative coexpression",
    y="sample",
    hue=hue,
    hue_order=unique_labels,
    palette=palette,
    cut=0,
    width=0.8,
    inner="quart",
)

plt.title(f"{plot_condition=} {plot_panel=} {method=} {target_count=}")
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), title=hue, frameon=False)
sns.despine(offset=10, trim=True)
plt.gca().xaxis.grid(True)
plt.savefig(out_file_plot_sample, dpi=300, bbox_inches="tight")
plt.close()


# Create panel violinplot
plt.figure(figsize=(6, 6))
g = sns.violinplot(
    data=df,
    x="log2 relative coexpression",
    y="panel",
    hue=hue,
    hue_order=unique_labels,
    palette=palette,
    cut=0,
    width=0.8,
)

plt.title(f"{plot_condition=} {method=} {target_count=}")
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5), title=hue, frameon=False)
sns.despine(offset=10, trim=True)
plt.gca().xaxis.grid(True)
plt.savefig(out_file_plot_panel, dpi=300, bbox_inches="tight")
plt.show()


df.to_parquet(out_file_gene_pairs)

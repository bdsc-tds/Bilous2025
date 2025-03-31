# %%
import dask

dask.config.set({"dataframe.query-planning": False})

import scanpy as sc
import numpy as np
import pandas as pd
import sys
import gc
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append("../../../workflow/scripts/")
import _utils
import readwrite

cfg = readwrite.config()
sns.set_style("ticks")


# %% [markdown]
# # Params

# %%
# cfg paths
xenium_dir = Path(cfg["xenium_processed_data_dir"])
xenium_count_correction_dir = Path(cfg["xenium_count_correction_dir"])
xenium_std_seurat_analysis_dir = Path(cfg["xenium_std_seurat_analysis_dir"])
xenium_cell_type_annotation_dir = Path(cfg["xenium_cell_type_annotation_dir"])
results_dir = Path(cfg["results_dir"])
palette_dir = Path(cfg["xenium_metadata_dir"])
std_seurat_analysis_dir = Path(cfg["xenium_std_seurat_analysis_dir"])
scrnaseq_processed_data_dir = Path(cfg["scrnaseq_processed_data_dir"])
seurat_to_h5_dir = results_dir / "seurat_to_h5"

# Params
signal_integrity_thresholds = [0.5, 0.7]
correction_methods = ["raw", "split_fully_purified", "resolvi", "resolvi_supervised"]
correction_methods += [
    f"ovrlpy_correction_{signal_integrity_threshold=}" for signal_integrity_threshold in signal_integrity_thresholds
]
num_samples = 30
mixture_k = 50
normalisation = "lognorm"
layer = "data"
layer_scrna = "RNA_counts"
reference = "matched_reference_combo"
method = "rctd_class_aware"
level = "Level2.1"
segmentation_palette = palette_dir / "col_palette_segmentation.csv"
count_correction_palette = palette_dir / "col_palette_correction_method.csv"

n_neighbors = 10
n_permutations = 30
n_repeats = 5
top_n = 20
scoring = "f1"
markers = "diffexpr"

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
plot_metrics = ["hypergeometric_pvalue", "NES", f"n_hits_{top_n=}", "mean_zscore_pvalue"]

labels_key = level

# %% [markdown]
# # Load corrected counts

# %%
xenium_paths = {}
xenium_annot_paths = {}

for correction_method in correction_methods:
    xenium_paths[correction_method] = {}
    xenium_annot_paths[correction_method] = {}

    for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
        if segmentation.stem == "proseg_mode":
            continue
        for condition in (conditions := segmentation.iterdir()):
            for panel in (panels := condition.iterdir()):
                for donor in (donors := panel.iterdir()):
                    for sample in (samples := donor.iterdir()):
                        k = (segmentation.stem, condition.stem, panel.stem, donor.stem, sample.stem)
                        name = "/".join(k)

                        # raw samples
                        if "proseg" in segmentation.stem:
                            k_proseg = ("proseg", condition.stem, panel.stem, donor.stem, sample.stem)
                            name_proseg = "/".join(k_proseg)
                            sample_dir = xenium_dir / f"{name_proseg}/raw_results"
                        else:
                            sample_dir = xenium_dir / f"{name}/normalised_results/outs"

                        sample_annotation = (
                            xenium_cell_type_annotation_dir
                            / f"{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet"
                        )

                        if correction_method == "raw":
                            xenium_paths[correction_method][k] = sample_dir
                            xenium_annot_paths[correction_method][k] = sample_annotation

                        # corrected samples
                        else:
                            if correction_method == "split_fully_purified":
                                name_corrected = f"{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/split_fully_purified/"
                                sample_corrected_counts_path = (
                                    xenium_count_correction_dir / f"{name_corrected}/corrected_counts.h5"
                                )

                            else:
                                if correction_method == "resolvi":
                                    name_corrected = f"{name}/{mixture_k=}/{num_samples=}/"
                                elif correction_method == "resolvi_supervised":
                                    name_corrected = f"{name}/{normalisation}/reference_based/{reference}/{method}/{level}/{mixture_k=}/{num_samples=}"
                                elif "ovrlpy" in correction_method:
                                    name_corrected = f"{name}"

                                sample_corrected_counts_path = (
                                    results_dir / f"{correction_method}/{name_corrected}/corrected_counts.h5"
                                )
                            sample_normalised_counts = (
                                xenium_std_seurat_analysis_dir
                                / f"{name}/{normalisation}/normalised_counts/{layer}.parquet"
                            )
                            sample_idx = (
                                xenium_std_seurat_analysis_dir
                                / f"{name}/{normalisation}/normalised_counts/cells.parquet"
                            )

                            xenium_paths[correction_method][k] = sample_corrected_counts_path


ads = readwrite.read_count_correction_samples(xenium_paths, correction_methods[1:])
ads["raw"] = readwrite.read_xenium_samples(xenium_paths["raw"], anndata=True, transcripts=False, max_workers=6)

# fix obs names for proseg expected, load cell types
# filter out cells without labels (this will apply QC thresholds as well since annotation is done after QC)
for correction_method in correction_methods:
    for k, ad in ads[correction_method].items():
        if ad is not None:
            if correction_method == "raw":
                if k[0] == "proseg_expected":
                    ad.obs_names = ad.obs_names.astype(str)
                    ad.obs_names = "proseg-" + ad.obs_names

                # filter cells and read labels for raw
                ad.obs[labels_key] = pd.read_parquet(xenium_annot_paths["raw"][k]).set_index("cell_id").iloc[:, 0]

                ad = ad[ad.obs[labels_key].notna()]
                if labels_key == "Level2.1":
                    # for custom Level2.1, simplify subtypes
                    ad.obs.loc[ad.obs[labels_key].str.contains("malignant"), labels_key] = "malignant cell"
                    ad.obs.loc[ad.obs[labels_key].str.contains("T cell"), labels_key] = "T cell"

                # remove tissue from cell type name
                ad.obs[labels_key] = ad.obs[labels_key].str.replace(r" of .+", "", regex=True)

                ads["raw"][k] = ad

            # filter cells and add labels from raw
            if correction_method != "raw":
                ad.obs[labels_key] = ads["raw"][k].obs[labels_key]
                ad = ad[[c for c in ads["raw"][k].obs_names if c in ad.obs_names]]
                ads[correction_method][k] = ad


# %%

# fix obs names for proseg expected, load cell types
# filter out cells without labels (this will apply QC thresholds as well since annotation is done after QC)
for correction_method in correction_methods:
    for k, ad in ads[correction_method].items():
        if ad is not None:
            if correction_method == "raw":
                if k[0] == "proseg_expected":
                    ad.obs_names = ad.obs_names.astype(str)
                    ad.obs_names = "proseg-" + ad.obs_names

                # filter cells and read labels for raw
                ad.obs[labels_key] = pd.read_parquet(xenium_annot_paths["raw"][k]).set_index("cell_id").iloc[:, 0]

                ad = ad[ad.obs[labels_key].notna()]
                if labels_key == "Level2.1":
                    # for custom Level2.1, simplify subtypes
                    ad.obs.loc[ad.obs[labels_key].str.contains("malignant"), labels_key] = "malignant cell"
                    ad.obs.loc[ad.obs[labels_key].str.contains("T cell"), labels_key] = "T cell"

                # remove tissue from cell type name
                ad.obs[labels_key] = ad.obs[labels_key].str.replace(r" of .+", "", regex=True)

                ads["raw"][k] = ad

            # filter cells and add labels from raw
            if correction_method != "raw":
                ad.obs[labels_key] = ads["raw"][k].obs[labels_key]
                ad = ad[[c for c in ads["raw"][k].obs_names if c in ad.obs_names]]
                ads[correction_method][k] = ad


# %% [markdown]
# # Load scRNA

# %%

ads_scrna = {}
for ref in (refs := scrnaseq_processed_data_dir.iterdir()):
    ref_name = ref.stem
    ref_dir = seurat_to_h5_dir / ref_name

    if "matched_combo_standard" not in ref_name:
        continue

    print(ref_name)

    ad = sc.read_10x_h5(ref_dir / f"{layer_scrna}.h5")
    ad.obs[labels_key] = pd.read_parquet(ref_dir / "metadata.parquet").set_index("cell_id")[labels_key]
    ad = ad[ad.obs[labels_key].notna()]

    if labels_key == "Level2.1":
        # for custom Level2.1, simplify subtypes
        ad.obs.loc[ad.obs[labels_key].str.contains("malignant"), labels_key] = "malignant cell"
        ad.obs.loc[ad.obs[labels_key].str.contains("T cell"), labels_key] = "T cell"
    # remove tissue from cell type name
    ad.obs[labels_key] = ad.obs[labels_key].str.replace(r" of .+", "", regex=True)

    ads_scrna[ref_name] = ad

ads_scrna["NSCLC"] = ads_scrna.pop("matched_combo_standard_breast_specific")
ads_scrna["breast"] = ads_scrna.pop("matched_combo_standard_lung_specific")


# Prepare pseudo bulk data
pbs_scrna = _utils.get_pseudo_bulk_scrna(ads_scrna, labels_key)
del ads_scrna
gc.collect()


# %% [markdown]
# # Plot cell type identity score from scRNA pseudobulk cosine similarity

# %%
cti = "B cell"
ref_panel = "lung"
plot_metric = "cosine_similarity"
palette = pd.read_csv(count_correction_palette, index_col=0).iloc[:, 0]

# get cell identity score df
df_all = _utils.get_cosine_similarity_score(
    ads,
    pbs_scrna,
    labels_key,
    correction_methods,
    columns=xenium_levels + ["correction_method", "cti", "cosine_similarity"],
)


_utils.rename_correction_methods(df_all)
df = df_all.query(f"panel == '{ref_panel}' and cti == '{cti}'")

# plotting params, palette
title = f"Cell type identity score for {cti=}"
unique_labels = [c for c in hue_correction_order if c in np.unique(df[hue_correction].dropna())]
unique_labels = unique_labels + [c for c in np.unique(df[hue_correction].dropna()) if c not in unique_labels]
palette = {u: palette[u] for u in unique_labels}
legend_handles = [mpatches.Patch(color=color, label=label) for label, color in palette.items()]

### hypergeometric pvalue boxplot
f = plt.figure(figsize=(10, 4))
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
)

sns.despine(offset=10, trim=True)
ax.yaxis.grid(True)

plt.suptitle(title)
f.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    title=hue_correction,
    frameon=False,
)
# plt.tight_layout(rect=[0, 0, 1, 0.95])
# plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
plt.show()

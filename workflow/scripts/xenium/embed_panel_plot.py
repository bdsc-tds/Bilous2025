from pathlib import Path
import sys
import pandas as pd
import scanpy as sc

sys.path.append("workflow/scripts/")
import readwrite

panel_path = Path(sys.argv[1])
embed_file = sys.argv[2]
out_file = sys.argv[3]

# read xenium samples
xenium_paths = {}
for sample in (samples := panel_path.iterdir()):
    for replicate in (replicates := sample.iterdir()):
        k = (sample.stem, replicate.stem)
        replicate_path = replicate / "normalised_results/outs"
        name = "/".join(k)

        xenium_paths[k] = replicate_path

# Read RCTD
rctd = {}
for k, path in xenium_paths.items():
    if (
        references := path.parents[1] / "cell_type_annotation/reference_based"
    ).exists():
        rctd[k] = {}

        for reference in (
            references := path.parents[1] / "cell_type_annotation/reference_based"
        ).iterdir():
            if reference.stem != "matched_reference":
                continue
            for method in (methods := reference.iterdir()):
                if method.stem != "rctd":
                    continue
                for level in (levels := method.iterdir()):
                    if level.stem != "Level2":
                        continue
                    cell_type_annotation_file = level / "single_cell/labels.csv"
                    if cell_type_annotation_file.exists():
                        rctd[k][reference.stem, method.stem, level.stem] = pd.read_csv(
                            cell_type_annotation_file, index_col=0
                        ).iloc[:, 0]


df = {}
for k in rctd:
    if len(rctd[k]):
        rctd_df = pd.DataFrame(rctd[k])
        rctd_df.columns = [col for col in rctd_df.columns]
        df[k] = rctd_df
df = pd.concat(df)  # .reset_index().dropna()
# df.columns = (*xenium_levels,'cell_id', *df.columns[6:])

obs = pd.read_parquet()

obs.index = pd.MultiIndex.from_tuples(
    obs["dataset_merge_id"].astype(object) + obs["cell_id"].map(lambda s: (s,))
)
obs.join(df)

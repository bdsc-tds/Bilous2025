import scipy
import pandas as pd
import argparse
from pathlib import Path


def transcripts_to_count_matrix(
    transcripts, cell_column="cell_id", feature_column="feature_name", qv_treshold=20
):
    transcripts = transcripts.query(
        f"(qv >= {qv_treshold}) & ({cell_column} != 'UNASSIGNED')"
    )
    cm = transcripts.pivot_table(
        index=cell_column, columns=feature_column, aggfunc="size", fill_value=0
    )
    return cm


# Set up argument parser
parser = argparse.ArgumentParser(description="Compute ovrlpy correction.")
parser.add_argument(
    "--sample_transcripts_path",
    type=Path,
    help="Path to the sample_signal_integrity file.",
)
parser.add_argument(
    "--sample_signal_integrity",
    type=Path,
    help="Path to the sample_signal_integrity file.",
)
parser.add_argument(
    "--out_file_corrected_counts",
    type=str,
    help="out_file_corrected_counts file to output.",
)
parser.add_argument(
    "--out_file_cells_mean_integrity",
    type=str,
    help="out_file_corrected_counts file to output.",
)
parser.add_argument(
    "--signal_integrity_threshold",
    type=float,
    help="signal_integrity_threshold parameter (threshold below which a pixel is low quality).",
)

args = parser.parse_args()

# Access the arguments
sample_transcripts_path = args.sample_transcripts_path
sample_signal_integrity = args.sample_signal_integrity
out_file_corrected_counts = args.out_file_corrected_counts
out_file_corrected_counts_index = args.out_file_corrected_counts_index
out_file_corrected_counts_columns = args.out_file_corrected_counts_columns
out_file_cells_mean_integrity = args.out_file_cells_mean_integrity
signal_integrity_threshold = args.signal_integrity_threshold


# load transcripts
coordinate_df = pd.read_parquet(sample_transcripts_path).query(
    "is_gene & (qv >= 20)"
)  # remove dummy & low qv molecules
coordinate_df["gene"] = coordinate_df["gene"].astype("category")


# load ovrlpy
signal_integrity = scipy.io.mmread(sample_signal_integrity).toarray()
# x and y are transposed in ovrlpy output. sanity check x and y must be transposed
assert (coordinate_df.y_pixel.max() + 1) == signal_integrity.shape[0]
assert (coordinate_df.x_pixel.max() + 1) == signal_integrity.shape[1]
coordinate_df["signal_integrity"] = signal_integrity[
    coordinate_df.y_pixel, coordinate_df.x_pixel
]

# filter transcripts based on ovrlpy signal integrity
coordinate_df_filtered = coordinate_df[
    coordinate_df.signal_integrity > signal_integrity_threshold
]

# compute mean cell integrity QC metric from unfiltered transcripts
cell_mean_integrity = (
    coordinate_df.query("cell_id != 'UNASSIGNED'")
    .groupby("cell_id")["signal_integrity"]
    .mean()
)
# cell_is_singlet = cell_mean_integrity > signal_integrity_threshold

# create filtered count matrix
corrected_counts = transcripts_to_count_matrix(
    coordinate_df_filtered, feature_column="gene"
).sum(1)


# store results
corrected_counts.to_parquet(out_file_corrected_counts)
cell_mean_integrity.to_frame().to_parquet(out_file_cells_mean_integrity)

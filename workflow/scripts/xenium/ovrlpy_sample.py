import ovrlpy
import scipy
import sys
import pandas as pd

replicate_transcripts_path = sys.argv[1]
out_file_signal_integrity = sys.argv[2]
out_file_signal_strength = sys.argv[3]
out_file_doublet_df = sys.argv[4]

# load data
coordinate_df = (
    pd.read_parquet(replicate_transcripts_path)
    .rename(
        columns={
            "x_location": "x",
            "y_location": "y",
            "z_location": "z",
            "feature_name": "gene",
        }
    )
    .query("is_gene & (qv >= 20)")
)  # remove dummy & low qv molecules

# run ovrlpy
signal_integrity, signal_strength, visualizer = ovrlpy.run(
    df=coordinate_df, cell_diameter=10, n_expected_celltypes=30
)
doublet_df = ovrlpy.detect_doublets(
    signal_integrity, signal_strength, minimum_signal_strength=3, integrity_sigma=2
)

# store results
signal_integrity = scipy.sparse.csr_matrix(signal_integrity)
signal_strength = scipy.sparse.csr_matrix(signal_strength)

scipy.io.mmwrite(out_file_signal_integrity, signal_integrity)
scipy.io.mmwrite(out_file_signal_strength, signal_strength)
doublet_df.to_parquet(out_file_doublet_df)

from pathlib import Path
import sys
import pandas as pd

sys.path.append("workflow/scripts/")
import coexpression
import readwrite

path = sys.argv[1]
out_file = sys.argv[2]
out_file_pos_rate = sys.argv[3]
method = sys.argv[4]

ad_ref = readwrite.read_xenium_sample(
    Path(path).stem, path, anndata_only=True, transcripts=False
)[1]
ad = readwrite.read_xenium_sample(
    Path(path).stem, path, anndata_only=True, transcripts=False
)[1]


CC = pd.read_parquet(out_file, index_col=0)
pos_rate = pd.read_parquet(out_file_pos_rate, index_col=0)["0"]

CCdiff, spurious_gene_pairs = coexpression.compare_segmentations(
    CCref,
    CCseg,
    pos_rate_ref_seg=pos_rate_ref_seg,
    pos_rate_other_seg=pos_rate_seg,
    min_positivity_rate=min_positivity_rate,
    cc_cutoff=cc_cutoff,
    method=method,
    log2=log2,
)

# CCdiff, spurious_gene_pairs

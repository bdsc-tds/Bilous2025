import pandas as pd
import scanpy as sc
import sys

sys.path.append("workflow/scripts/")
import readwrite
import coexpression

# params
sample_counts = sys.argv[1]
out_file = sys.argv[2]
out_file_pos_rate = sys.argv[3]
method = sys.argv[4]
target_count = int(sys.argv[5])
gene_panel_path = sys.argv[6]

# read counts
adata = sc.read_10x_h5(sample_counts)

# read gene panel names
gene_panel = readwrite.get_gene_panel_info(gene_panel_path)["name"]

# subset adata
adata = adata[:, [gene for gene in gene_panel if gene in adata.var_names]]
if adata.shape[1] == 0:
    raise ValueError("No genes from gene panel found in adata.var_names")

# compute coexpression
CC, X_downsampled, pos, pos_rate, mask = coexpression.coexpression(adata, target_count=target_count, method=method)

# save as parquet
CC.to_parquet(out_file)
pos_rate.to_frame().to_parquet(out_file_pos_rate)

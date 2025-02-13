import scanpy as sc
import sys

sys.path.append("workflow/scripts/")
import coexpression

# params
path = sys.argv[1]
out_file = sys.argv[2]
out_file_pos_rate = sys.argv[3]
method = sys.argv[4]
target_count = int(sys.argv[5])

# read counts
adata = sc.read_10x_h5(path)

# compute coexpression
CC, X_downsampled, pos, pos_rate, mask = coexpression.coexpression(adata, target_count=target_count, method=method)

# save as parquet
CC.to_parquet(out_file)
pos_rate.to_frame().to_parquet(out_file_pos_rate)

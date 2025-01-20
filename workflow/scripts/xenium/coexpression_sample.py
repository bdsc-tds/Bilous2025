from pathlib import Path
import sys
sys.path.append('workflow/scripts/')
import coexpression
import readwrite

path = sys.argv[1]
out_file = sys.argv[2]
out_file_pos_rate = sys.argv[3]
method = sys.argv[4]
target_count = int(sys.argv[5])

ad = readwrite.read_xenium_sample(Path(path).stem,path,anndata_only=True,transcripts=False)[1]
CC, X_downsampled, pos, pos_rate, mask = coexpression.coexpression(ad,target_count=target_count,method=method)

CC.to_parquet(out_file)
pos_rate.to_frame().to_parquet(out_file_pos_rate)
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


# cfg paths
xenium_dir = Path(cfg['xenium_processed_data_dir'])
results_dir = Path(cfg['results_dir'])

# Segmentation, mapping paths
dir_segmentations = {
    dir_segmentation.name: (dir_segmentation)
    for dir_segmentation in xenium_dir.iterdir()
}

# Read resegmentations and RCTD

xenium_paths = {}
umaps = {}
cc_paths = []

for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 
        for panel in (panels := cohort.iterdir()):
            for sample in (samples := panel.iterdir()):
                for replicate in (replicates := sample.iterdir()):
                    
                    k = (segmentation.stem,cohort.stem,panel.stem,sample.stem,replicate.stem)
                    replicate_path = replicate / "normalised_results/outs"
                    name = '/'.join(k)

                    xenium_paths[k] = replicate_path


xenium_levels = ('segmentation','cohort','panel','sample','replicate')
ads = readwrite.read_xenium_samples(xenium_paths,anndata_only=True,transcripts=False,sample_name_as_key=False)
ads = pd.Series(ads.values(),
                index=pd.Index(ads.keys(),name = xenium_levels),
                dtype=object).sort_index()


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

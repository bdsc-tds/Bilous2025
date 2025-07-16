import pandas as pd

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
scrnaseq_processed_data_dir = Path(config['scrnaseq_processed_data_dir'])
seurat_to_h5_dir = results_dir / 'seurat_to_h5'

# params from pipeline config
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# Params
layer = 'RNA_counts'

n_comps = config['umap_n_comps']
n_neighbors = config['umap_n_neighbors']
min_dist = config['umap_min_dist']
metric = config['umap_metric']

reference1 = scrnaseq_processed_data_dir / 'matched_combo_standard_lung_specific'
reference1_name = reference1.stem
reference1_dir = seurat_to_h5_dir / reference1_name
reference1_is_done = reference1_dir / '.done'

reference2 = scrnaseq_processed_data_dir / 'matched_combo_standard_breast_specific'
reference2_name = reference2.stem
reference2_dir = seurat_to_h5_dir / reference2_name
reference2_is_done = reference2_dir / '.done'

out_file = results_dir / f'embed_panel_scrnaseq_joint/{reference1_name}_{reference2_name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 

rule embed_panel_scrnaseq_joint:
    input:
        reference1_is_done=reference1_is_done,
        reference2_is_done=reference2_is_done
    output:
        out_file=out_file,
    params:
        reference1=reference1_dir,
        reference2=reference2_dir,
        layer=layer,
        n_comps=n_comps,
        n_neighbors=n_neighbors,
        metric=metric,
        min_dist=min_dist,
        min_counts=min_counts,
        min_features=min_features,
        max_counts=max_counts,
        max_features=max_features,
        min_cells=min_cells,
    threads: 1
    resources:
        mem='50GB',
        runtime='30m',
        # slurm_partition = "gpu",
        # slurm_extra = '--gres=gpu:1',
    conda:
        "general_cuda"
    shell:
        """
        mkdir -p "$(dirname {output.out_file})"

        python -u workflow/scripts/scRNAseq/embed_panel_scrnaseq_joint.py \
            --reference1 {params.reference1} \
            --reference2 {params.reference2} \
            --out_file {output.out_file} \
            --layer {params.layer} \
            --n_comps {params.n_comps} \
            --n_neighbors {params.n_neighbors} \
            --metric {params.metric} \
            --min_dist {params.min_dist} \
            --min_counts {params.min_counts} \
            --min_features {params.min_features} \
            --max_counts {params.max_counts} \
            --max_features {params.max_features} \
            --min_cells {params.min_cells} \
            
        echo "DONE"
        """
# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# params from pipeline config
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# Params
n_comps = 50
n_neighbors = 50
min_dist = 0.3
metric = 'cosine'

out_files_panel = []
for segmentation in (segmentations := xenium_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):

            k = (segmentation.stem,condition.stem,panel.stem)
            name = '/'.join(k)

            out_file = results_dir / f'embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
            out_files_panel.append(out_file)

            rule:
                name: f'embed_panel/{name}'
                input:
                    panel=panel,
                output:
                    out_file=out_file,
                params:
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
                    mem='100GB' if panel.stem == '5k' else '50GB',
                    runtime='30m' if panel.stem == '5k' else '20m',
                    slurm_partition = "gpu",
                    slurm_extra = '--gres=gpu:1',
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file})"

                    python workflow/scripts/xenium/embed_panel.py \
                        --panel {input.panel} \
                        --out_file {output.out_file} \
                        --n_comps {params.n_comps} \
                        --n_neighbors {params.n_neighbors} \
                        --metric {params.metric} \
                        --min_dist {params.min_dist} \
                        --min_counts {params.min_counts} \
                        --min_features {params.min_features} \
                        --max_counts {params.max_counts} \
                        --max_features {params.max_features} \
                        --min_cells {params.min_cells}

                    echo "DONE"
                    """


out_files_condition = []
for segmentation in (segmentations := xenium_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 

            k = (segmentation.stem,condition.stem)
            name = '/'.join(k)

            out_file = results_dir / f'embed_condition/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
            out_files_condition.append(out_file)

            rule:
                name: f'embed_condition/{name}'
                input:
                    condition=condition,
                output:
                    out_file=out_file,
                params:
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
                    mem='100GB' if panel.stem == '5k' else '50GB',
                    runtime='30m' if panel.stem == '5k' else '20m',
                    slurm_partition = "gpu",
                    slurm_extra = '--gres=gpu:1'
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file})"

                    python workflow/scripts/xenium/embed_condition.py \
                        --condition {input.condition} \
                        --out_file {output.out_file} \
                        --n_comps {params.n_comps} \
                        --n_neighbors {params.n_neighbors} \
                        --metric {params.metric} \
                        --min_dist {params.min_dist} \
                        --min_counts {params.min_counts} \
                        --min_features {params.min_features} \
                        --max_counts {params.max_counts} \
                        --max_features {params.max_features} \
                        --min_cells {params.min_cells}

                    echo "DONE"
                    """

rule embed_panel_all:
    input:
        out_files_panel

rule embed_condition_all:
    input:
        out_files_condition
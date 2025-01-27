# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
n_comps = 50
n_neighbors = 30
min_dist = 0.3
metric = 'cosine'

out_files = []

for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 
        for panel in (panels := cohort.iterdir()):

            k = (segmentation.stem,cohort.stem,panel.stem)
            name = '/'.join(k)

            if replicate_transcripts_path.exists():

                out_file = results_dir / f'embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric=}.parquet' 
                out_files.append(out_file)

            rule:
                name: f'embed_panel/{name}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric=}'
                input:
                    panel_path=panel.stem,
                output:
                    out_file=out_file,
                params:
                    n_comps=n_comps,
                    metric=metric,
                    min_dist=min_dist,
                threads: 1
                resources:
                    mem='100GB' if panel.stem == '5k' else '50GB',
                    runtime='1h' if panel.stem == '5k' else '30m',
                    partition = "gpu"
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file})"

                    python workflow/scripts/scRNAseq/embed_panel.py \
                    {input.panel_path} \
                    {output.out_file} \
                    {params.n_comps}
                    {params.metric} \
                    {params.min_dist} \

                    echo "DONE"
                    """



rule embed_panel_samples:
    input:
        out_files
# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])

# Params
references = ['matched_reference','external_reference']
methods = ['rctd']
levels = ['Level2']

out_files = []

for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 
        for panel in (panels := cohort.iterdir()):

            k = (segmentation.stem,cohort.stem,panel.stem)
            name = '/'.join(k)
            embed_file = results_dir / f'embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric=}.parquet'

            for reference in references:
                for method in methods:
                    for level in levels:

                        out_file = figures_dir / f"embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric=}_{reference=}_{method=}_{level=}.png"
                        out_files.append(out_file)

                    rule:
                        name: f'embed_panel_plot/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric=}_{reference=}_{method=}_{level=}'
                        input:
                            panel_path=panel.stem,
                            embed_file=embed_file,
                        output:
                            out_file=out_file,
                        params:
                            references=references
                        threads: 1
                        resources:
                            mem='30GB',
                            runtime='10m',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file})"

                            python workflow/scripts/scRNAseq/embed_panel_plot.py \
                            {input.panel_path} \
                            {input.embed_file} \
                            {output.out_file} \

                            echo "DONE"
                            """

rule embed_panel_samples_plot:
    input:
        out_files
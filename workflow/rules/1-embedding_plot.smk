# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])

# Params
n_comps = 50
n_neighbors = 150
min_dist = 0.3
metric = 'cosine'

references = ['matched_reference','external_reference']
methods = ['rctd']
levels = ['Level2','cohort','replicate',] # cohort and replicate as color to plot added here in addition to levels
extension = 'png'

out_files_panel = []

for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 
        for panel in (panels := cohort.iterdir()):

            k = (segmentation.stem,cohort.stem,panel.stem)
            name = '/'.join(k)
            embed_file = results_dir / f'embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet'

            for reference in references:
                for method in methods:
                    for level in levels:
                        if level == 'cohort':
                            continue

                        out_file = figures_dir / f"embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}_{reference}_{method}_{level}.{extension}"
                        out_files_panel.append(out_file)

                        rule:
                            name: f'embed_panel_plot/{name}/umap_{reference}_{method}_{level}'
                            input:
                                panel=panel,
                                embed_file=embed_file,
                            output:
                                out_file=out_file,
                            params:
                                reference=reference,
                                method=method,
                                level=level,
                            threads: 1
                            resources:
                                mem='30GB',
                                runtime='10m',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file})"

                                python workflow/scripts/xenium/embed_panel_plot.py \
                                --panel {input.panel} \
                                --embed_file {input.embed_file} \
                                --reference {params.reference} \
                                --method {params.method} \
                                --level {params.level} \
                                --out_file {output.out_file} \

                                echo "DONE"
                                """



out_files_cohort = []
for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 
        for panel in (panels := cohort.iterdir()):

            k = (segmentation.stem,cohort.stem)
            name = '/'.join(k)
            embed_file = results_dir / f'embed_cohort/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet'

            for reference in references:
                for method in methods:
                    for level in levels:

                        out_file = figures_dir / f"embed_cohort/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}_{reference}_{method}_{level}.{extension}"
                        out_files_cohort.append(out_file)

                        rule:
                            name: f'embed_cohort_plot/{name}/umap_{reference}_{method}_{level}'
                            input:
                                cohort=cohort,
                                embed_file=embed_file,
                            output:
                                out_file=out_file,
                            params:
                                reference=reference,
                                method=method,
                                level=level,
                            threads: 1
                            resources:
                                mem='30GB',
                                runtime='10m',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file})"

                                python workflow/scripts/xenium/embed_cohort_plot.py \
                                --cohort {input.cohort} \
                                --embed_file {input.embed_file} \
                                --reference {params.reference} \
                                --method {params.method} \
                                --level {params.level} \
                                --out_file {output.out_file} \

                                echo "DONE"
                                """





rule embed_panel_samples_plot:
    input:
        out_files_panel

rule embed_cohort_samples_plot:
    input:
        out_files_cohort
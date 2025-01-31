from pathlib import Path

# cfg paths
silhouette_dir = Path(config['results_dir']) / 'silhouette'
figures_dir = Path(config['figures_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
segmentation_palette = palette_dir / 'col_palette_segmentation.csv'
references = ['matched_reference','external_reference']
methods = ['rctd_class_aware']
levels = ['Level2'] 
extension = 'png'

out_files = []
for segmentation in (segmentations := silhouette_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for reference in references:
                for method in methods:
                    for level in levels:

                        k = (segmentation.stem,condition.stem,panel.stem)
                        name = '/'.join(k)

                        out_file = figures_dir / f'silhouette_plot/{name}/silhouette_{reference}_{method}_{level}.{extension}'
                        out_files.append(out_file)

                        rule:
                            name: f'silhouette_plot/{name}/silhouette_{reference}_{method}_{level}'
                            input:
                                panel=panel,
                            output:
                                out_file=out_file,
                            params:
                                segmentation_palette=segmentation_palette,
                                reference=reference,
                                method=method,
                                level=level,
                            threads: 1
                            resources:
                                mem='5GB',
                                runtime='5m',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file})"

                                python workflow/scripts/xenium/silhouette_panel_plot.py \
                                --panel {input.panel} \
                                --out_file {output.out_file} \
                                --segmentation_palette {params.segmentation_palette} \
                                --reference {params.reference} \
                                --method {params.method} \
                                --level {params.level} \

                                echo "DONE"
                                """


rule silhouette_plot_donors:
    input:
        out_files
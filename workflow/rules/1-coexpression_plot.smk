from pathlib import Path

# cfg paths
silhouette_dir = Path(config['results_dir']) / 'silhouette'
figures_dir = Path(config['figures_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
methods = ['conditional','jaccard','pearson','spearman']
target_counts = [30,50,200]
ref_segmentation = '10x_0um'
ref_oversegmentation = '10x_15um'
segmentation_palette = palette_dir / 'col_palette_segmentation.csv'
extension = 'png'

out_files = []
for segmentation in (segmentations := silhouette_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for method in methods:
                for target_count in target_counts:
                    if target_count > 50 and panel.stem != '5k':
                        continue

                        k = (segmentation.stem,condition.stem,panel.stem)
                        name = '/'.join(k)

                        out_file = figures_dir / f'coexpression_plot/{name}/coexpression_{method=}_{target_count=}.{extension}'
                        out_files.append(out_file)

                        rule:
                            name: f'coexpression_plot_panel/{name}/coexpression_{method=}_{target_count=}'
                            input:
                                panel=panel,
                            output:
                                out_file=out_file,
                            params:
                                reference=reference,
                                method=method,
                                target_count=target_count,
                                min_positivity_rate=min_positivity_rate,
                                cc_cutoff=cc_cutoff,
                                log2=log2,
                                ref_segmentation=ref_segmentation,
                                ref_oversegmentation=ref_oversegmentation,
                                segmentation_palette=segmentation_palette,
                            threads: 1
                            resources:
                                mem='5GB',
                                runtime='5m',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file})"

                                python workflow/scripts/xenium/coexpression_panel_plot.py \
                                --panel {input.panel} \
                                --out_file {output.out_file} \
                                --reference {params.reference} \
                                --method {params.method} \
                                --target_count {params.target_count} \
                                --min_positivity_rate {params.min_positivity_rate} \
                                --cc_cutoff {params.cc_cutoff} \
                                --log2 {params.log2} \
                                --ref_segmentation {params.ref_segmentation} \
                                --ref_oversegmentation {params.ref_oversegmentation} \
                                --segmentation_palette {params.segmentation_palette} \

                                echo "DONE"
                                """


rule coexpression_plot_panel:
    input:
        out_files




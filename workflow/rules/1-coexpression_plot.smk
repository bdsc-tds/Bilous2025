from pathlib import Path

# cfg paths
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
methods = ['conditional','jaccard','pearson','spearman']
target_counts = [30,50,200]
min_positivity_rate = 0.01
cc_cutoff = 1.5
log2 = True
ref_segmentation = '10x_0um'
ref_oversegmentation = '10x_15um'
segmentation_palette = palette_dir / 'col_palette_segmentation.csv'
extension = 'png'

out_files = []
for segmentation in (segmentations := (results_dir/'coexpression').iterdir()):
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

                    out_file_plot = figures_dir / f'coexpression_plot/{name}/coexpression_{method}_{target_count=}.{extension}'
                    out_file_gene_pairs = results_dir / f'coexpression_gene_pairs/{name}/coexpression_gene_pairs_{method}_{target_count=}.parquet'
                    out_files.extend([out_file_plot,out_file_gene_pairs])

                    rule:
                        name: f'coexpression_plot_panel/{name}/coexpression_{method}_{target_count=}'
                        input:
                            coexpression_is_done=results_dir / "coexpression.done",
                            panel=panel,
                        output:
                            out_file_plot=out_file_plot,
                            out_file_gene_pairs=out_file_gene_pairs,
                        params:
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
                            mem='80GB' if panel.stem == '5k' else '20GB',
                            runtime='40m' if panel.stem == '5k' else '10m',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file_plot})"

                            python workflow/scripts/xenium/coexpression_panel_plot.py \
                            --panel {input.panel} \
                            --out_file_plot {output.out_file_plot} \
                            --out_file_gene_pairs {output.out_file_gene_pairs} \
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


rule coexpression_plot_all:
    input:
        out_files




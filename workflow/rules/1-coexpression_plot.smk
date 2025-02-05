from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
methods = ['conditional','jaccard']#,'pearson','spearman']
target_counts = [30,50,200]
min_positivity_rate = 0.01
# min_cond_coex = 0.05
cc_cutoff = 1.5
ref_segmentation = '10x_0um'
ref_oversegmentation = '10x_15um'
segmentation_palette = palette_dir / 'col_palette_segmentation.csv'
extension = 'png'

out_files = []
for segmentation in (segmentations := xenium_dir.iterdir()):
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

                    panel_coexpression = results_dir / f'coexpression/{name}'
                    out_file_plot_sample = figures_dir / f'coexpression_plot/{name}/coexpression_{method}_{target_count=}_sample.{extension}'
                    out_file_plot_panel = figures_dir / f'coexpression_plot/{name}/coexpression_{method}_{target_count=}_panel.{extension}'
                    out_file_gene_pairs = results_dir / f'coexpression_gene_pairs/{name}/coexpression_gene_pairs_{method}_{target_count=}.parquet'
                    out_files.extend([out_file_plot_sample,out_file_plot_panel,out_file_gene_pairs])

                    # adapt resources
                    if panel.stem == '5k':
                        if target_count > 50:
                            mem = '120GB'
                            runtime = '40m'
                        else:
                            mem = '80GB'
                            runtime = '40m'
                    else:
                        mem = '20GB'
                        runtime = '10m'

                    rule:
                        name: f'coexpression_plot_panel/{name}/coexpression_{method}_{target_count=}'
                        input:
                            coexpression_is_done=results_dir / "coexpression.done",
                        output:
                            out_file_plot_sample=out_file_plot_sample,
                            out_file_plot_panel=out_file_plot_panel,
                            out_file_gene_pairs=out_file_gene_pairs,
                        params:
                            panel=panel_coexpression,
                            method=method,
                            target_count=target_count,
                            min_positivity_rate=min_positivity_rate,
                            # min_cond_coex=min_cond_coex if method in ['conditional','jaccard'] else 0.,
                            cc_cutoff=cc_cutoff,
                            log2=True if method in ['conditional','jaccard'] else False,
                            ref_segmentation=ref_segmentation,
                            ref_oversegmentation=ref_oversegmentation,
                            segmentation_palette=segmentation_palette,
                        threads: 1
                        resources:
                            mem=mem,
                            runtime=runtime,
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file_plot_sample})"

                            python workflow/scripts/xenium/coexpression_panel_plot.py \
                            --panel {params.panel} \
                            --out_file_plot_sample {output.out_file_plot_sample} \
                            --out_file_plot_panel {output.out_file_plot_panel} \
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
                            # --min_cond_coex {params.min_cond_coex} \


rule coexpression_plot_all:
    input:
        out_files




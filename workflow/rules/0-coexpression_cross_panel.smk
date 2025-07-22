# Params
genes = 'conditions'
methods = ['conditional','jaccard']#,'pearson','spearman']
target_counts = [15]
out_files = []

for segmentation in (segmentations := xenium_processed_data_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):

                    if segmentation.stem == 'proseg':
                        sample_path = sample / 'raw_results'
                        segmentation_name = 'proseg_expected'
                    else:
                        sample_path = sample / "normalised_results/outs"
                        segmentation_name = segmentation.stem
                    
                    k = (segmentation_name,condition.stem,panel.stem,donor.stem,sample.stem)
                    name = '/'.join(k)

                    if sample_path.exists():

                        for method in methods:
                            for target_count in target_counts:

                                out_file_coexpr = results_dir / f'coexpression_cross_panel/{name}/coexpression_{method}_{target_count}.parquet' 
                                out_file_pos_rate = results_dir / f'coexpression_cross_panel/{name}/positivity_rate_{method}_{target_count}.parquet'

                                out_files.extend([out_file_coexpr,out_file_pos_rate])

                                rule:
                                    name: f'coexpression_cross_panel/{name}/{method}_{target_count}'
                                    input:
                                        sample_path=sample_path,
                                    output:
                                        out_file_coexpr=out_file_coexpr,
                                        out_file_pos_rate=out_file_pos_rate,
                                    params:
                                        method=method,
                                        target_count=target_count,
                                        genes=genes,
                                    threads: 1
                                    resources:
                                        mem='40GB' if panel.stem == '5k' else '20GB',
                                        runtime='20m' if panel.stem == '5k' else '10m',
                                    conda:
                                        "spatial"
                                    shell:
                                        """
                                        mkdir -p "$(dirname {output.out_file_coexpr})"

                                        python workflow/scripts/coexpression_sample.py \
                                        -i {input.sample_path} \
                                        -m {params.method} \
                                        -c {params.target_count} \
                                        --outcoexp {output.out_file_coexpr} \
                                        --outposrate {output.out_file_pos_rate} \
                                        -g {params.genes}

                                        echo "DONE"
                                        """


rule coexpression_cross_panel_all:
    input:
        out_files
    output:
        touch(results_dir / "coexpression_cross_panel.done")

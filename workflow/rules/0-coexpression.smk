from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
methods = ['conditional','jaccard']#,'pearson','spearman']
target_counts = [30,50,200]
out_files = []

for segmentation in (segmentations := xenium_dir.iterdir()):
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
                                if target_count > 50 and panel.stem != '5k':
                                    continue

                                out_file_coexpr = results_dir / f'coexpression/{name}/coexpression_{method}_{target_count}.parquet' 
                                out_file_pos_rate = results_dir / f'coexpression/{name}/positivity_rate_{method}_{target_count}.parquet'

                                out_files.extend([out_file_coexpr,out_file_pos_rate])

                                rule:
                                    name: f'coexpression/{name}/{method}_{target_count}'
                                    input:
                                        sample_path=sample_path,
                                    output:
                                        out_file_coexpr=out_file_coexpr,
                                        out_file_pos_rate=out_file_pos_rate,
                                    params:
                                        method=method,
                                        target_count=target_count,
                                    threads: 1
                                    resources:
                                        mem='40GB' if panel.stem == '5k' else '20GB',
                                        runtime='20m' if panel.stem == '5k' else '10m',
                                    conda:
                                        "spatial"
                                    shell:
                                        """
                                        mkdir -p "$(dirname {output.out_file_coexpr})"

                                        python workflow/scripts/xenium/coexpression_sample.py \
                                        {input.sample_path} \
                                        {output.out_file_coexpr} \
                                        {output.out_file_pos_rate} \
                                        {params.method} \
                                        {params.target_count} \

                                        echo "DONE"
                                        """


rule coexpression_all:
    input:
        out_files
    output:
        touch(results_dir / "coexpression.done")

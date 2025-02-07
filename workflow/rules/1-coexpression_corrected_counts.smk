from pathlib import Path

# cfg paths
# xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_std_seurat_analysis_dir =  Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])

# Params
methods = ['conditional','jaccard','pearson','spearman']
target_counts = [30,50,200]
out_files_resolvi = []

for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):

                    k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                    name = '/'.join(k)
                    sample_path_resolvi = results_dir / f'resolvi/{name}/resolvi_corrected.parquet'

                    if sample_path.exists():

                        for method in methods:
                            for target_count in target_counts:
                                if target_count > 50 and panel.stem != '5k':
                                    continue

                                out_file_coexpr = results_dir / f'resolvi_coexpression/{name}/coexpression_{method}_{target_count}.parquet' 
                                out_file_pos_rate = results_dir / f'resolvi_coexpression/{name}/positivity_rate_{method}_{target_count}.parquet'

                                out_files_resolvi.extend([out_file_coexpr,out_file_pos_rate])

                                rule:
                                    name: f'resolvi_coexpression/{name}/{method}_{target_count}'
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

                                        python workflow/scripts/xenium/coexpression_sample_resolvi.py \
                                        {input.sample_path} \
                                        {output.out_file_coexpr} \
                                        {output.out_file_pos_rate} \
                                        {params.method} \
                                        {params.target_count} \

                                        echo "DONE"
                                        """


rule coexpression_resolvi_all:
    input:
        out_files_resolvi
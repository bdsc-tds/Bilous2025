from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
signal_integrity_thresholds = [0.5,0.7]
correction_methods = ['resolvi'] + [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
methods = ['conditional','jaccard']#,'pearson','spearman']
target_counts = [30,50,200]
out_files = []

for correction_method in correction_methods:
    for segmentation in (segmentations := xenium_dir.iterdir()):
        if segmentation.stem == 'proseg_v1':
            continue
        for condition in (conditions := segmentation.iterdir()): 
            for panel in (panels := condition.iterdir()):
                for donor in (donors := panel.iterdir()):
                    for sample in (samples := donor.iterdir()):
                        if donor.stem in ['0WMU','1G73']:
                            continue
                        k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                        name = '/'.join(k)

                        sample_path = results_dir / f'{correction_method}/{name}/corrected_counts.h5'

                    # if sample_path.exists():

                        for method in methods:
                            for target_count in target_counts:
                                if target_count > 50 and panel.stem != '5k':
                                    continue

                                out_file_coexpr = results_dir / f'{correction_method}_coexpression/{name}/coexpression_{method}_{target_count}.parquet' 
                                out_file_pos_rate = results_dir / f'{correction_method}_coexpression/{name}/positivity_rate_{method}_{target_count}.parquet'

                                out_files.extend([out_file_coexpr,out_file_pos_rate])

                                rule:
                                    name: f'{correction_method}_coexpression/{name}/{method}_{target_count}'
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
                                        mem='60GB' if panel.stem == '5k' else '20GB',
                                        runtime='20m' if panel.stem == '5k' else '10m',
                                    conda:
                                        "spatial"
                                    shell:
                                        """
                                        mkdir -p "$(dirname {output.out_file_coexpr})"

                                        python workflow/scripts/coexpression_sample.py \
                                        -i {input.sample_path} \
                                        --outcoexp {output.out_file_coexpr} \
                                        --outposrate {output.out_file_pos_rate} \
                                        -m {params.method} \
                                        -c {params.target_count} \

                                        echo "DONE"
                                        """


rule coexpression_corrected_counts_all:
    input:
        out_files
    output:
        touch([results_dir / f"{correction_method}_coexpression.done"
            for correction_method in correction_methods]
        )
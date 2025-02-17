from pathlib import Path

# cfg paths
# xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
correction_methods = ['resolvi','ovrlpy_correction']
methods = ['conditional','jaccard']#,'pearson','spearman']
target_counts = [30,50,200]
signal_integrity_threshold = 0.5
out_files = []

for correction_method in correction_methods:
    for segmentation in (segmentations := (results_dir / correction_method).iterdir()):
        for condition in (conditions := segmentation.iterdir()): 
            for panel in (panels := condition.iterdir()):
                for donor in (donors := panel.iterdir()):
                    for sample in (samples := donor.iterdir()):

                        k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                        name = '/'.join(k)

                        if correction_method == 'resolvi':
                            sample_path = results_dir / f'{correction_method}/{name}/resolvi_corrected_counts.h5'
                        elif correction_method == 'ovrlpy_correction':
                            sample_path = results_dir / f'{correction_method}/{name}/corrected_counts_{signal_integrity_threshold=}.h5'
                        else:
                            raise ValueError(f"{correction_method} not supported. Use 'resolvi' or 'ovrlpy_correction'")

                        if sample_path.exists():

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
                                            mem='40GB' if panel.stem == '5k' else '20GB',
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
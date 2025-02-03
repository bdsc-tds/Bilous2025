from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
signal_integrity_threshold = 0.5

segmentation = xenium_dir / '10x_0um' # ovrlpy output does not depend on segmentation, just run for 10x_0um
out_files = []

for condition in (conditions := segmentation.iterdir()): 
    for panel in (panels := condition.iterdir()):
        for donor in (donors := panel.iterdir()):
            for sample in (samples := donor.iterdir()):

                k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                sample_transcripts_path = sample / "normalised_results/outs/transcripts.parquet"
                name = '/'.join(k)

                # if sample_transcripts_path.exists():

                out_file_signal_integrity = results_dir / f'ovrlpy/{name}/signal_integrity.mm' 
                out_file_signal_strength = results_dir / f'ovrlpy/{name}/signal_strength.mm'
                out_file_doublet_df = results_dir / f'ovrlpy/{name}/doublet_df.parquet'

                out_files.extend([out_file_signal_integrity,
                                    out_file_signal_strength,
                                    out_file_doublet_df])

                rule:
                    name: f'ovrlpy/{name}'
                    input:
                        sample_transcripts_path=sample_transcripts_path,
                    output:
                        out_file_signal_integrity=out_file_signal_integrity,
                        out_file_signal_strength=out_file_signal_strength,
                        out_file_doublet_df=out_file_doublet_df,
                    threads: 1
                    resources:
                        mem='500GB' if panel.stem == '5k' else '200GB',
                        runtime='10h' if panel.stem == '5k' else '1h',
                    conda:
                        "spatial"
                    shell:
                        """
                        mkdir -p "$(dirname {output.out_file_signal_integrity})"

                        python workflow/scripts/xenium/ovrlpy_sample.py \
                        --sample_transcripts_path {input.sample_transcripts_path} \
                        --out_file_signal_integrity {output.out_file_signal_integrity} \
                        --out_file_signal_strength {output.out_file_signal_strength} \
                        --out_file_doublet_df {output.out_file_doublet_df} \

                        echo "DONE"
                        """

for segmentation in (segmentations := xenium_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):

                    k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                    name = '/'.join(k)
                    sample_transcripts_path = sample / "normalised_results/outs/transcripts.parquet"
                    sample_signal_integrity = results_dir / f'ovrlpy/{name}/signal_integrity.mm'

                    # if sample_transcripts_path.exists():
                    out_file_corrected_counts = results_dir / f'ovrlpy_correction/{name}/corrected_counts_{signal_integrity_threshold=}.mm'
                    out_file_corrected_counts_index = results_dir / f'ovrlpy_correction/{name}/corrected_counts_{signal_integrity_threshold=}_index.mm'
                    out_file_corrected_counts_columns = results_dir / f'ovrlpy_correction/{name}/corrected_counts_{signal_integrity_threshold=}_columns.mm'
                    out_file_cells_mean_integrity = results_dir / f'ovrlpy_correction/{name}/cells_mean_integrity.parquet'

                    out_files.extend([
                        out_file_corrected_counts,
                        out_file_corrected_counts_index,
                        out_file_corrected_counts_columns,
                        out_file_cells_mean_integrity,])

                    rule:
                        name: f'ovrlpy/{name}'
                        input:
                            sample_transcripts_path=sample_transcripts_path,
                            sample_signal_integrity=sample_signal_integrity,
                        output:
                            out_file_corrected_counts=out_file_corrected_counts,
                            out_file_corrected_counts_index=out_file_corrected_counts_index,
                            out_file_corrected_counts_columns=out_file_corrected_counts_columns,
                            out_file_cells_mean_integrity=out_file_cells_mean_integrity,
                        params:
                            signal_integrity_threshold=signal_integrity_threshold
                        threads: 1
                        resources:
                            mem='30GB' if panel.stem == '5k' else '20GB',
                            runtime='20m',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file_corrected_counts})"

                            python workflow/scripts/xenium/ovrlpy_sample.py \
                            --sample_transcripts_path {input.sample_transcripts_path} \
                            --sample_signal_integrity {input.sample_signal_integrity} \
                            --out_file_corrected_counts {output.out_file_corrected_counts} \
                            --out_file_corrected_counts_index {output.out_file_corrected_counts_index} \
                            --out_file_corrected_counts_columns {output.out_file_corrected_counts_columns} \
                            --out_file_cells_mean_integrity {output.out_file_cells_mean_integrity} \
                            --signal_integrity_threshold {params.signal_integrity_threshold} \

                            echo "DONE"
                            """



rule ovrlpy_all:
    input:
        out_files
from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_raw_data_dir = Path(config['xenium_raw_data_dir'])
results_dir = Path(config['results_dir'])

# Params
segmentation = xenium_dir / '10x_0um' # ovrlpy output does not depend on segmentation, just run for 10x_0um
out_files = []

for cohort in (cohorts := segmentation.iterdir()): 
    for panel in (panels := cohort.iterdir()):
        for sample in (samples := panel.iterdir()):
            for replicate in (replicates := sample.iterdir()):

                k = (segmentation.stem,cohort.stem,panel.stem,sample.stem,replicate.stem)
                replicate_transcripts_path = replicate / "normalised_results/outs/transcripts.parquet"
                name = '/'.join(k)

                if replicate_transcripts_path.exists():

                    out_file_signal_integrity = results_dir / f'ovrlpy/{name}/signal_integrity.mm' 
                    out_file_signal_strength = results_dir / f'ovrlpy/{name}/signal_strength.mm'
                    out_file_doublet_df = results_dir / f'ovrlpy/{name}/doublet_df.parquet'

                    out_files.extend([out_file_signal_integrity,
                                        out_file_signal_strength,
                                        out_file_doublet_df])

                    rule:
                        name: f'ovrlpy/{name}'
                        input:
                            replicate_transcripts_path=replicate_transcripts_path,
                        output:
                            out_file_signal_integrity=out_file_signal_integrity,
                            out_file_signal_strength=out_file_signal_strength,
                            out_file_doublet_df=out_file_doublet_df,
                        threads: 1
                        resources:
                            mem='300GB' if panel.stem == '5k' else '200GB',
                            runtime='3h' if panel.stem == '5k' else '1h',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file_signal_integrity})"

                            python workflow/scripts/xenium/ovrlpy_sample.py \
                            {input.replicate_transcripts_path} \
                            {output.out_file_signal_integrity} \
                            {output.out_file_signal_strength} \
                            {output.out_file_doublet_df} \

                            echo "DONE"
                            """


rule ovrlpy_samples:
    input:
        out_files
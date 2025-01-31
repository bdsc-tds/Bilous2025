from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
segmentation = xenium_dir / '10x_0um' # ovrlpy output does not depend on segmentation, just run for 10x_0um
out_files = []

for condition in (conditions := segmentation.iterdir()): 
    for panel in (panels := condition.iterdir()):
        for donor in (donors := panel.iterdir()):
            for sample in (samples := donor.iterdir()):

                k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                sample_transcripts_path = sample / "normalised_results/outs/transcripts.parquet"
                name = '/'.join(k)

                if sample_transcripts_path.exists():

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
                            mem='400GB' if panel.stem == '5k' else '200GB',
                            runtime='10h' if panel.stem == '5k' else '1h',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file_signal_integrity})"

                            python workflow/scripts/xenium/ovrlpy_donor.py \
                            {input.sample_transcripts_path} \
                            {output.out_file_signal_integrity} \
                            {output.out_file_signal_strength} \
                            {output.out_file_doublet_df} \

                            echo "DONE"
                            """


rule ovrlpy_donors:
    input:
        out_files
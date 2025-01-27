from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
max_sample_size = 50_000
out_files = []

for segmentation in (segmentations := xenium_dir.iterdir()):
    for cohort in (cohorts := segmentation.iterdir()): 
        for panel in (panels := cohort.iterdir()):
            for sample in (samples := panel.iterdir()):
                for replicate in (replicates := sample.iterdir()):

                    k = (segmentation.stem,cohort.stem,panel.stem,sample.stem,replicate.stem)
                    replicate_path = replicate / "normalised_results/outs"
                    name = '/'.join(k)

                    if replicate_path.exists():

                        out_file = results_dir / f'silhouette/{name}/doublet_df.parquet'
                        out_files.append(out_file)

                        rule:
                            name: f'silhouette/{name}'
                            input:
                                replicate_path=replicate_path,
                            output:
                                out_file=out_file,
                            params:
                                max_sample_size=max_sample_size
                            threads: 1
                            resources:
                                mem='100GB',
                                runtime='1h',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file})"

                                python workflow/scripts/xenium/silhouette_sample.py \
                                {input.replicate_path} \
                                {output.out_file} \
                                {params.max_sample_size}

                                echo "DONE"
                                """


rule silhouette_samples:
    input:
        out_files
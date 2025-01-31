from pathlib import Path
import yaml

f = '/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/containers/skang/xenium_analysis_pipeline/samples/small_samples.yml'
with open(f) as f:
    donors_config = yaml.safe_load(f)

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])

# Params
max_donor_size = 1_000
out_files = []

for segmentation in (segmentations := xenium_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    if panel.stem not in donors_config[condition.stem]:
                        continue
                    if donor.stem not in donors_config[condition.stem][panel.stem]:
                        continue

                    k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                    sample_path = sample / "normalised_results/outs"
                    name = '/'.join(k)

                    if sample_path.exists():

                        out_file = results_dir / f'silhouette/{name}/silhouette.parquet'
                        out_files.append(out_file)

                        rule:
                            name: f'silhouette/{name}'
                            input:
                                sample_path=sample_path,
                            output:
                                out_file=out_file,
                            params:
                                max_donor_size=max_donor_size
                            threads: 1
                            resources:
                                mem='30GB',
                                runtime='15m',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file})"

                                python workflow/scripts/xenium/silhouette_donor.py \
                                {input.sample_path} \
                                {output.out_file} \
                                {params.max_donor_size}

                                echo "DONE"
                                """


rule silhouette_donors:
    input:
        out_files
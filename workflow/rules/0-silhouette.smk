from pathlib import Path
import yaml

f = '/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/containers/skang/xenium_analysis_pipeline/samples/small_samples.yml'
with open(f) as f:
    donors_config = yaml.safe_load(f)

# cfg paths
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
xenium_cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
results_dir = Path(config['results_dir'])

# Params
normalisations = ['lognorm','sctransform']
max_sample_size = 50_000
out_files = []

for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    if panel.stem not in donors_config[condition.stem]:
                        continue
                    if donor.stem not in donors_config[condition.stem][panel.stem]:
                        continue

                    for normalisation in normalisations:

                        k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,normalisation)
                        name = '/'.join(k)

                        sample_counts = sample / f'{normalisation}/normalised_counts/counts.parquet'
                        sample_idx = sample / f'{normalisation}/normalised_counts/cells.parquet'
                        sample_annotation_dir = xenium_cell_type_annotation_dir / f'{name}/reference_based'

                        # sample_path = sample / "normalised_results/outs"


                        if sample_path.exists():

                            out_file = results_dir / f'silhouette/{name}/silhouette.parquet'
                            out_files.append(out_file)

                            rule:
                                name: f'silhouette/{name}'
                                input:
                                    sample_counts=sample_counts,
                                    sample_idx=sample_idx,
                                    sample_annotation_dir=sample_annotation_dir,
                                output:
                                    out_file=out_file,
                                params:
                                    max_sample_size=max_sample_size
                                threads: 1
                                resources:
                                    mem='30GB',
                                    runtime='15m',
                                conda:
                                    "spatial"
                                shell:
                                    """
                                    mkdir -p "$(dirname {output.out_file})"

                                    python workflow/scripts/xenium/silhouette_sample.py \
                                    {input.sample_path} \
                                    {input.sample_idx} \
                                    {input.sample_annotation_dir} \
                                    {output.out_file} \
                                    {params.max_sample_size}

                                    echo "DONE"
                                    """


rule silhouette_all:
    input:
        out_files
    output:
        touch(results_dir / "silhouette.done")
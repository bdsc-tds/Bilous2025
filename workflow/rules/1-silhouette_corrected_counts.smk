from pathlib import Path
import yaml

# f = '/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/containers/skang/xenium_analysis_pipeline/samples/small_samples.yml'
# with open(f) as f:
#     donors_config = yaml.safe_load(f)

# cfg paths
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
xenium_cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
results_dir = Path(config['results_dir'])

# Params
metric = 'euclidean'
normalisations = ['lognorm','sctransform']
layers = ['data','scale_data']
max_sample_size = 50_000
out_files = []

for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    # if panel.stem not in donors_config[condition.stem]:
                    #     continue
                    # if donor.stem not in donors_config[condition.stem][panel.stem]:
                    #     continue

                    for normalisation in normalisations:
                        for layer in layers:

                            k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,normalisation)
                            name = '/'.join(k)

                            sample_pca = sample / f'{normalisation}/preprocessed/pca.parquet'
                            sample_idx = sample / f'{normalisation}/preprocessed/cells.parquet'
                            sample_annotation_dir = xenium_cell_type_annotation_dir / f'{name}/reference_based'

                            out_file = results_dir / f'silhouette/{name}/silhouette_{layer}.parquet'
                            out_files.append(out_file)

                            rule:
                                name: f'silhouette/{name}_{layer}'
                                input:
                                    sample_pca=sample_pca,
                                    sample_idx=sample_idx,
                                    sample_annotation_dir=sample_annotation_dir,
                                output:
                                    out_file=out_file,
                                params:
                                    max_sample_size=max_sample_size,
                                    metric=metric,
                                threads: 1
                                resources:
                                    mem='30GB',
                                    runtime='15m',
                                conda:
                                    "general_cuda"
                                shell:
                                    """
                                    mkdir -p "$(dirname {output.out_file})"

                                    python workflow/scripts/xenium/silhouette_sample.py \
                                    --sample_pca {input.sample_pca} \
                                    --sample_idx {input.sample_idx} \
                                    --sample_annotation_dir {input.sample_annotation_dir} \
                                    --out_file {output.out_file} \
                                    --max_sample_size {params.max_sample_size} \
                                    --metric {params.metric}

                                    echo "DONE"
                                    """


rule silhouette_all:
    input:
        out_files
    output:
        touch(results_dir / "silhouette.done")
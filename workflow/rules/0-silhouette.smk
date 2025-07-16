from pathlib import Path
import yaml

# cfg paths
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
xenium_cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
results_dir = Path(config['results_dir'])

# Params
normalisations = ['lognorm','sctransform']
layers = ['data','scale_data']
max_sample_size = 50_000
metric = 'euclidean'

out_files = []
for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):

                    for normalisation in normalisations:
                        for layer in layers:

                            k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,normalisation)
                            name = '/'.join(k)
                            name_annot = '/'.join(k[:-1]+('lognorm',))


                            sample_pca = sample / f'{normalisation}/preprocessed/pca.parquet'
                            sample_idx = sample / f'{normalisation}/preprocessed/cells.parquet'
                            sample_annotation_dir = xenium_cell_type_annotation_dir / f'{name_annot}/reference_based'

                            if not sample_annotation_dir.exists():
                                continue
                                
                            out_file = results_dir / f'silhouette/{name}/silhouette_{layer}_{metric}.parquet'
                            out_files.append(out_file)

                            rule:
                                name: f'silhouette/{name}_{layer}_{metric}'
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
                                    runtime='30m',
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
                                    --metric {params.metric} \

                                    echo "DONE"
                                    """


rule silhouette_all:
    input:
        out_files
    output:
        touch(results_dir / "silhouette.done")
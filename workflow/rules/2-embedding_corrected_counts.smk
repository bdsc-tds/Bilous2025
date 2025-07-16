# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
xenium_count_correction_dir = Path(config['xenium_count_correction_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
results_dir = Path(config['results_dir'])

# params from pipeline config
signal_integrity_thresholds = [0.5,0.7]
# correction_methods = ['split_fully_purified','resolvi','resolvi_supervised'] + [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
correction_methods = ['split_fully_purified']
normalisations = ['lognorm']
layers = ['data']#,'scale_data']
references = ['matched_reference_combo','external_reference']
methods = ['rctd_class_aware']
levels = ['Level2.1']
cell_type_normalisation = 'lognorm'

min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# Params
n_comps = config['umap_n_comps']
n_neighbors = config['umap_n_neighbors']
min_dist = config['umap_min_dist']
metric = config['umap_metric']
raw_corrected_counts = True

# resolvi params
num_samples = 30
mixture_k = 50

out_files_panel = []

for correction_method in correction_methods:
    for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
        if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
            continue
        for condition in (conditions := segmentation.iterdir()): 
            for panel in (panels := condition.iterdir()):
                for reference in references:
                    if reference == 'external_reference' and condition.stem == 'mesothelioma_pilot':
                        continue
                    for method in methods:
                        for level in levels:
                            if level == 'Level2.1' and reference == 'external_reference':
                                continue
                            for normalisation in normalisations: 
                                for layer in layers:
                                    k = (segmentation.stem,condition.stem,panel.stem)
                                    name = '/'.join(k)

                                    if correction_method == 'split_fully_purified':
                                        panel_path = xenium_count_correction_dir / name
                                    else:
                                        panel_path = results_dir / f'{correction_method}/{name}'     
                                                                                                
                                    out_file = results_dir / f'{correction_method}_embed_panel/{name}/{normalisation}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                                    out_files_panel.append(out_file)

                                    rule:
                                        name: f'{correction_method}_embed_panel/{name}/{normalisation}_{layer}_{reference}'
                                        input:
                                            count_correction_is_done=results_dir / f"{correction_method}.done"
                                        output:
                                            out_file=out_file,
                                        params:
                                            panel=panel_path,
                                            normalisation=normalisation,
                                            layer=layer,
                                            reference=reference,
                                            method=method,
                                            level=level,
                                            n_comps=n_comps,
                                            n_neighbors=n_neighbors,
                                            metric=metric,
                                            min_dist=min_dist,
                                            min_counts=min_counts,
                                            min_features=min_features,
                                            max_counts=max_counts,
                                            max_features=max_features,
                                            min_cells=min_cells,
                                            num_samples=num_samples,
                                            mixture_k=mixture_k,
                                            xenium_count_correction_dir=xenium_count_correction_dir,
                                            results_dir=results_dir,
                                            correction_method=correction_method,
                                            cell_type_annotation_dir = cell_type_annotation_dir,
                                            cell_type_normalisation = cell_type_normalisation,
                                            raw_corrected_counts='--raw_corrected_counts' if raw_corrected_counts else '',
                                        threads: 1
                                        resources:
                                            mem='100GB' if panel.stem == '5k' else '50GB',
                                            # runtime='30m' if panel.stem == '5k' else '20m',
                                            runtime='8h',
                                            # slurm_partition = "gpu",
                                            # slurm_extra = '--gres=gpu:1',
                                        conda:
                                            "general_cuda"
                                        shell:
                                            """
                                            mkdir -p "$(dirname {output.out_file})"

                                            python workflow/scripts/xenium/embed_panel_corrected_counts.py \
                                                --panel {params.panel} \
                                                --out_file {output.out_file} \
                                                --normalisation {params.normalisation} \
                                                --layer {params.layer} \
                                                --n_comps {params.n_comps} \
                                                --n_neighbors {params.n_neighbors} \
                                                --metric {params.metric} \
                                                --min_dist {params.min_dist} \
                                                --min_counts {params.min_counts} \
                                                --min_features {params.min_features} \
                                                --max_counts {params.max_counts} \
                                                --max_features {params.max_features} \
                                                --min_cells {params.min_cells} \
                                                --num_samples {params.num_samples} \
                                                --mixture_k {params.mixture_k} \
                                                --xenium_count_correction_dir {params.xenium_count_correction_dir} \
                                                --results_dir {params.results_dir} \
                                                --correction_method {params.correction_method} \
                                                --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                                                --cell_type_normalisation {params.cell_type_normalisation} \
                                                --reference {params.reference} \
                                                --method {params.method} \
                                                --level {params.level} \
                                                {params.raw_corrected_counts}
                                                
                                            echo "DONE"
                                            """

rule embed_panel_corrected_counts_all:
    input:
        out_files_panel
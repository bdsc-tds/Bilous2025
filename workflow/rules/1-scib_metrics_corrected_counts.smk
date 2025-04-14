# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
xenium_count_correction_dir = Path(config['xenium_count_correction_dir'])

# Params
n_comps = 50
max_n_cells = 100_000

normalisations = ['lognorm']#,'sctransform']
layers = ['data']#,'scale_data']
references = ['matched_reference_combo']#,'external_reference']
methods = ['rctd_class_aware']
levels = ['Level2.1'] # condition and sample as color to plot added here in addition to levels


# params from pipeline config
signal_integrity_thresholds = [0.5,0.7]
correction_methods = ['split_fully_purified','resolvi','resolvi_supervised'] + [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
normalisations = ['lognorm']
layers = ['data','scale_data']
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# Params
raw_corrected_counts = True

# resolvi params
num_samples = 30
mixture_k = 50


out_files_panel = []

for correction_method in correction_methods:
    for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
        if segmentation.stem in ['proseg_mode']:
            continue
        # if correction_method == 'split_fully_purified' and segmentation.stem in ['10x_mm_5um']:
        #     continue
        for condition in (conditions := segmentation.iterdir()): 
            for panel in (panels := condition.iterdir()):
                for normalisation in normalisations:
                    for layer in layers: 
                        for reference in references:
                            for method in methods:
                                for level in levels:

                                    # input embedding file (doesn't depend on ref,method or color loops but more readable to have here)
                                    k = (segmentation.stem,condition.stem,panel.stem)
                                    name = '/'.join(k)
                                    
                                    if correction_method == 'split_fully_purified':
                                        panel_path = xenium_count_correction_dir / name
                                    else:
                                        panel_path = results_dir / f'{correction_method}/{name}'

                                    out_file = results_dir / f"scib_metrics_panel/{correction_method}/{name}/{normalisation}/scib_metrics_{layer}_{reference}_{method}_{level}_{n_comps=}_{max_n_cells=}.parquet"
                                    out_files_panel.append(out_file)

                                    rule:
                                        name: f'scib_metrics_panel/{correction_method}/{name}/scib_metrics_{layer}_{reference}_{method}_{level}'
                                        input:
                                            panel=panel_path,
                                        output:
                                            out_file=out_file,
                                        params:
                                            cell_type_annotation_dir=cell_type_annotation_dir,
                                            normalisation=normalisation,
                                            layer=layer,
                                            reference=reference,
                                            method=method,
                                            level=level,
                                            n_comps=n_comps,
                                            min_counts=min_counts,
                                            min_features=min_features,
                                            max_counts=max_counts,
                                            max_features=max_features,
                                            min_cells=min_cells,
                                            max_n_cells=max_n_cells,
                                            num_samples=num_samples,
                                            mixture_k=mixture_k,
                                            xenium_count_correction_dir=xenium_count_correction_dir,
                                            results_dir=results_dir,
                                            correction_method=correction_method,
                                            raw_corrected_counts='--raw_corrected_counts' if raw_corrected_counts else '',
                                        threads: 1
                                        resources:
                                            mem='80GB',
                                            runtime='5h',
                                        conda:
                                            "spatial"
                                        shell:
                                            """
                                            mkdir -p "$(dirname {output.out_file})"

                                            python workflow/scripts/xenium/scib_metrics_panel_corrected_counts.py \
                                            --panel {input.panel} \
                                            --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                                            --out_file {output.out_file} \
                                            --normalisation {params.normalisation} \
                                            --layer {params.layer} \
                                            --reference {params.reference} \
                                            --method {params.method} \
                                            --level {params.level} \
                                            --n_comps {params.n_comps} \
                                            --max_n_cells {params.max_n_cells} \
                                            --max_counts {params.max_counts} \
                                            --max_features {params.max_features} \
                                            --min_cells {params.min_cells} \
                                            --num_samples {params.num_samples} \
                                            --mixture_k {params.mixture_k} \
                                            --xenium_count_correction_dir {params.xenium_count_correction_dir} \
                                            --results_dir {params.results_dir} \
                                            --correction_method {params.correction_method} \
                                            {params.raw_corrected_counts}

                                            echo "DONE"
                                            """



rule scib_metrics_panel_corrected_counts_all:
    input:
        out_files_panel
    output:
        touch(results_dir / "scib_metrics_corrected_counts.done")
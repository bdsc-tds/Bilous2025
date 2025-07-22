segmentation = sorted(xenium_processed_data_dir.iterdir())[0] # arbitrary segmentation just to loop over conditions and panels

params_product = list(product(normalisations, layers, references, methods, levels))

out_files = []

for condition in (conditions := segmentation.iterdir()): 
    for panel in (panels := condition.iterdir()):
        for normalisation, layer, reference, method, level in params_product:
            if level not in CONDITIONS_LEVELS[condition.stem]:
                continue
            if reference not in CONDITIONS_REFERENCES[condition.stem]:
                continue

            k = (condition.stem,panel.stem)
            name = '/'.join(k)

            out_dir = figures_dir / f'contamination_metrics_cosine_similarity_scrnaseq_boxplot/{name}/{normalisation}/{layer}_{reference}_{method}_{level}/'
            out_file = out_dir / '.done'
            out_files.append(out_file)

            rule:
                name: f'contamination_metrics_cosine_similarity_scrnaseq_boxplot/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                input:
                    resolvi_is_done = results_dir / "resolvi.done",
                    resolvi_supervised_is_done = results_dir / "resolvi_supervised.done",
                    ovrlpy_is_done = [results_dir / f"ovrlpy_correction_{signal_integrity_threshold=}.done"  for signal_integrity_threshold in signal_integrity_thresholds]
                output:
                    out_file = touch(out_file)
                params:
                    condition=condition.stem,
                    panel=panel.stem,
                    correction_methods=correction_methods,
                    results_dir=results_dir,
                    xenium_processed_data_dir=xenium_processed_data_dir,
                    xenium_count_correction_dir=xenium_count_correction_dir,
                    scrnaseq_processed_data_dir=scrnaseq_processed_data_dir,
                    seurat_to_h5_dir=seurat_to_h5_dir,
                    std_seurat_analysis_dir=std_seurat_analysis_dir,
                    cell_type_annotation_dir=cell_type_annotation_dir,
                    out_dir=out_dir,
                    normalisation=normalisation,
                    layer=layer,
                    reference=reference,
                    method=method,
                    level=level,
                    mixture_k=mixture_k,
                    num_samples=num_samples,
                    use_precomputed="--use_precomputed" if use_precomputed else "",
                    count_correction_palette=count_correction_palette,
                    dpi=dpi,
                    extension=extension,
                threads: 1
                resources:
                    mem='100GB',
                    runtime='30m',
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file})"

                    python workflow/scripts/xenium/contamination_metrics_diffexpr_cosine_similarity_boxplot.py \
                    --condition {params.condition} \
                    --panel {params.panel} \
                    --correction_methods {params.correction_methods} \
                    --results_dir {params.results_dir} \
                    --xenium_processed_data_dir {params.xenium_processed_data_dir} \
                    --xenium_count_correction_dir {params.xenium_count_correction_dir} \
                    --scrnaseq_processed_data_dir {params.scrnaseq_processed_data_dir} \
                    --seurat_to_h5_dir {params.seurat_to_h5_dir} \
                    --std_seurat_analysis_dir {params.std_seurat_analysis_dir} \
                    --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                    --out_dir {params.out_dir} \
                    --normalisation {params.normalisation} \
                    --layer {params.layer} \
                    --reference {params.reference} \
                    --method {params.method} \
                    --level {params.level} \
                    --mixture_k {params.mixture_k} \
                    --num_samples {params.num_samples} \
                    --use_precomputed {params.use_precomputed} \
                    --count_correction_palette {params.count_correction_palette} \
                    --dpi {params.dpi} \
                    --extension {params.extension} \

                    echo "DONE"
                    """


rule contamination_metrics_cosine_similarity_scrnaseq_boxplot_all:
    input:
        out_files


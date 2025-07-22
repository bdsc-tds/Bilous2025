segmentation = sorted(xenium_processed_data_dir.iterdir())[0] # arbitrary segmentation just to loop over conditions and panels

params_product = list(product(normalisations, layers, references, methods, levels))

out_files = []
for markers_mode in markers_modes:

    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for normalisation, layer, reference, method, level in params_product:
                if level not in CONDITIONS_LEVELS[condition.stem]:
                    continue
                if reference not in CONDITIONS_REFERENCES[condition.stem]:
                    continue

                k = (condition.stem,panel.stem)
                name = '/'.join(k)
                name_params = f"{markers_mode}_{radius=}_{top_n=}"

                out_dir = figures_dir / f'contamination_metrics_{name_params}_specificity_boxplot/{name}/{normalisation}/{layer}_{reference}_{method}_{level}/'
                out_file = out_dir / '.done'
                out_files.append(out_file)

                rule:
                    name: f'contamination_metrics_{name_params}_specificity_boxplot/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                    input:
                        contamination_metrics_is_done=results_dir / f"contamination_metrics_{name_params}.done",
                        contamination_metrics_corrected_counts_is_done=results_dir / f"contamination_metrics_{name_params}_corrected_counts.done",
                    output:
                        out_file = touch(out_file)
                    params:
                        condition=condition.stem,
                        panel=panel.stem,
                        correction_methods=correction_methods,
                        results_dir=results_dir,
                        std_seurat_analysis_dir=std_seurat_analysis_dir,
                        cell_type_annotation_dir=cell_type_annotation_dir,
                        out_dir=out_dir,
                        normalisation=normalisation,
                        layer=layer,
                        reference=reference,
                        method=method,
                        level=level,
                        top_n=top_n,
                        mixture_k=mixture_k,
                        num_samples=num_samples,
                        use_precomputed="--use_precomputed" if use_precomputed else "",
                        count_correction_palette=count_correction_palette,
                        radius=radius,
                        dpi=dpi,
                        extension=extension,
                    threads: 1
                    resources:
                        mem='30GB',
                        runtime='2h',
                    conda:
                        "spatial"
                    shell:
                        """
                        mkdir -p "$(dirname {output.out_file})"

                        python workflow/scripts/xenium/contamination_metrics_diffexpr_specificity_boxplot.py \
                            --condition {params.condition} \
                            --panel {params.panel} \
                            --correction_methods {params.correction_methods} \
                            --results_dir {params.results_dir} \
                            --std_seurat_analysis_dir {params.std_seurat_analysis_dir} \
                            --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                            --out_dir {params.out_dir} \
                            --normalisation {params.normalisation} \
                            --layer {params.layer} \
                            --reference {params.reference} \
                            --method {params.method} \
                            --level {params.level} \
                            --top_n {params.top_n} \
                            --mixture_k {params.mixture_k} \
                            --num_samples {params.num_samples} \
                            {params.use_precomputed} \
                            --count_correction_palette {params.count_correction_palette} \
                            --radius {params.radius} \
                            --dpi {params.dpi} \
                            --extension {params.extension} \

                        echo "DONE"
                        """


rule contamination_metrics_specificity_boxplot_all:
    input:
        out_files


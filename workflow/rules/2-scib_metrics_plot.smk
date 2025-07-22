biocons_scores = [
    "cLISI",
    "Isolated labels",
    "KMeans NMI",
    "KMeans ARI",
    "Leiden NMI",
    "Leiden ARI",
    "Silhouette label",
    "Bio conservation",
]
batchcor_scores = ["iLISI", "Graph connectivity", "KBET", "Silhouette batch", "Batch correction"]
scores = biocons_scores + batchcor_scores

segmentation = sorted(xenium_processed_data_dir.iterdir())[0] # arbitrary segmentation just to loop over conditions and panels

params_product = list(product(normalisations, layers, references, methods, levels, scores))

out_files_panel = []

for condition in (conditions := segmentation.iterdir()): 
    for panel in (panels := condition.iterdir()):
        for normalisation, layer, reference, method, level, score in params_product:

            if level not in CONDITIONS_LEVELS[condition.stem]:
                continue
            if reference not in CONDITIONS_REFERENCES[condition.stem]:
                continue

            # input embedding file (doesn't depend on ref,method or color loops but more readable to have here)
            k = (condition.stem,panel.stem,normalisation)
            name = '/'.join(k)

            name_score = score.replace(" ","_")
            out_file = figures_dir / f"scib_metrics_panel_plot/{name}/scib_metrics_{name_score}_{layer}_{reference}_{method}_{level}_{n_comps=}_{max_n_cells=}.{extension}"
            out_files_panel.append(out_file)

            rule:
                name: f'scib_metrics_panel_plot/{name}/scib_metrics_{name_score}_{layer}_{reference}_{method}_{level}'
                input:
                    scib_metrics_is_done=results_dir / "scib_metrics.done",
                    scib_metrics_corrected_counts_is_done=results_dir / "scib_metrics_corrected_counts.done",
                output:
                    out_file=out_file,
                params:
                    condition=condition.stem,
                    panel=panel.stem,
                    correction_methods=correction_methods,
                    std_seurat_analysis_dir=std_seurat_analysis_dir,
                    cell_type_annotation_dir=cell_type_annotation_dir,
                    scib_metrics_results_dir=scib_metrics_results_dir,
                    out_file=out_file,
                    normalisation=normalisation,
                    layer=layer,
                    reference=reference,
                    method=method,
                    level=level,
                    n_comps=n_comps,
                    max_n_cells=max_n_cells,
                    count_correction_palette=count_correction_palette,
                    dpi=dpi,
                    score=score,
                threads: 1
                resources:
                    mem='80GB',
                    runtime='5h',
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file})"

                    python workflow/scripts/xenium/scib_metrics_panel_plot.py \
                    --condition {params.condition} \
                    --panel {params.panel} \
                    --correction_methods {params.correction_methods} \
                    --std_seurat_analysis_dir {params.std_seurat_analysis_dir} \
                    --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                    --scib_metrics_results_dir {params.scib_metrics_results_dir} \
                    --out_file {params.out_file} \
                    --normalisation {params.normalisation} \
                    --layer {params.layer} \
                    --reference {params.reference} \
                    --method {params.method} \
                    --level {params.level} \
                    --n_comps {params.n_comps} \
                    --max_n_cells {params.max_n_cells} \
                    --count_correction_palette {params.count_correction_palette} \
                    --dpi {params.dpi} \
                    --score "{params.score}" \

                    echo "DONE"
                    """



rule scib_metrics_panel_plot_all:
    input:
        out_files_panel

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
xenium_count_correction_dir = Path(config['xenium_count_correction_dir'])
scib_metrics_results_dir = results_dir / "scib_metrics_panel"
palette_dir = Path(config['xenium_metadata_dir'])

# Params
normalisations = ['lognorm']#,'sctransform']
layers = ['data']#,'scale_data']
references = ['matched_reference_combo','external_reference']
methods = ['rctd_class_aware']
levels = ['Level2.1'] # condition and sample as color to plot added here in addition to levels
dpi = 300
extension = 'png'
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

# params from pipeline config
n_comps = 50
max_n_cells = 100_000
signal_integrity_thresholds = [0.5,0.7]
correction_methods = ['raw','split_fully_purified','resolvi','resolvi_supervised'] + [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5
count_correction_palette = palette_dir / 'col_palette_correction_method.csv'

# Params
raw_corrected_counts = True

# resolvi params
num_samples = 30
mixture_k = 50
segmentation = sorted(xenium_dir.iterdir())[0] # arbitrary segmentation just to loop over conditions and panels


out_files_panel = []

for condition in (conditions := segmentation.iterdir()): 
    for panel in (panels := condition.iterdir()):
        for normalisation in normalisations:
            for layer in layers: 
                for reference in references:
                    for method in methods:
                        for level in levels:
                            if level == 'Level2.1' and reference == 'external_reference':
                                continue
                            if level in ['Level1','Level2'] and reference == 'matched_reference_combo':
                                continue
                            for score in scores:

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
                                        "general_cuda"
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

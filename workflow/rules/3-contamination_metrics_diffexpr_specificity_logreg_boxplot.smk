from pathlib import Path
import yaml
import itertools
import pandas as pd

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])
palette_dir = Path(config['xenium_metadata_dir'])
figures_dir = Path(config['figures_dir'])

# Params
# probably only need to run for lognorm data
signal_integrity_thresholds = [0.5,0.7]
correction_methods = ['raw','split_fully_purified','resolvi','resolvi_supervised'] + [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
normalisations = ['lognorm',]
layers = ['data',]
references = ['matched_reference_combo']
methods = ['rctd_class_aware']
levels = ['Level2.1']
use_precomputed = True

radius = 15
n_permutations = 30
n_splits= 5
n_repeats = 5
top_n = 20
scoring = 'precision'
cv_mode = 'spatial'
markers_modes = ['diffexpr']#,'common_markers'] #'/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/markers/cellmarker_cell_types_markers.json'

# resolvi params
num_samples = 30
mixture_k = 50

dpi = 300
extension = 'png'

segmentation = sorted(xenium_dir.iterdir())[0] # arbitrary segmentation just to loop over conditions and panels
count_correction_palette = palette_dir / 'col_palette_correction_method.csv'

out_files = []
for markers_mode in markers_modes:

    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for normalisation in normalisations:
                for layer in layers:
                    for reference in references:
                        for method in methods:
                            for level in levels:

                                k = (condition.stem,panel.stem)
                                name = '/'.join(k)
                                name_params = f"{markers_mode}_{radius=}_{n_permutations=}_{n_splits=}_{top_n=}_{scoring}_{cv_mode}"

                                out_dir = figures_dir / f'contamination_metrics_{name_params}_specificity_logreg_boxplot/{name}/{normalisation}/{layer}_{reference}_{method}_{level}/'
                                out_file = out_dir / '.done'
                                out_files.append(out_file)

                                rule:
                                    name: f'contamination_metrics_{name_params}_specificity_logreg_boxplot/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                                    input:
                                        contamination_metrics_is_done=results_dir / f"contamination_metrics_{name_params}_logreg.done",
                                        contamination_metrics_corrected_counts_is_done=results_dir / f"contamination_metrics_{name_params}_logreg_corrected_counts.done",
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
                                        cv_mode=cv_mode,
                                        n_splits=n_splits,
                                        n_permutations=n_permutations,
                                        n_repeats=n_repeats,
                                        scoring=scoring,
                                        dpi=dpi,
                                        extension=extension,
                                    threads: 1
                                    resources:
                                        mem='30GB',
                                        runtime='30m',
                                    conda:
                                        "spatial"
                                    shell:
                                        """
                                        mkdir -p "$(dirname {output.out_file})"

                                        python workflow/scripts/xenium/contamination_metrics_diffexpr_specificity_logreg_boxplot.py \
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
                                            --cv_mode {params.cv_mode} \
                                            --n_splits {params.n_splits} \
                                            --n_permutations {params.n_permutations} \
                                            --n_repeats {params.n_repeats} \
                                            --scoring {params.scoring} \
                                            --dpi {params.dpi} \
                                            --extension {params.extension} \

                                        echo "DONE"
                                        """


rule contamination_metrics_specificity_logreg_boxplot_all:
    input:
        out_files


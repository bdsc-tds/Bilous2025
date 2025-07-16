from pathlib import Path
import yaml
import itertools
import pandas as pd

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
count_correction_dir = Path(config['xenium_count_correction_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
scrnaseq_processed_data_dir = Path(config['scrnaseq_processed_data_dir'])
results_dir = Path(config['results_dir'])
palette_dir = Path(config['xenium_metadata_dir'])
figures_dir = Path(config['figures_dir'])
seurat_to_h5_dir = results_dir / 'seurat_to_h5'

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
dpi = 300
extension = 'png'

# resolvi params
num_samples = 30
mixture_k = 50

segmentation = sorted(xenium_dir.iterdir())[0] # arbitrary segmentation just to loop over conditions and panels
count_correction_palette = palette_dir / 'col_palette_correction_method.csv'

out_files = []

for condition in (conditions := segmentation.iterdir()): 
    for panel in (panels := condition.iterdir()):
        for normalisation in normalisations:
            for layer in layers:
                for reference in references:
                    for method in methods:
                        for level in levels:

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
                                    xenium_dir=xenium_dir,
                                    count_correction_dir=count_correction_dir,
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
                                    "general_cuda"
                                shell:
                                    """
                                    mkdir -p "$(dirname {output.out_file})"

                                    python workflow/scripts/xenium/contamination_metrics_diffexpr_cosine_similarity_boxplot.py \
                                    --condition {params.condition} \
                                    --panel {params.panel} \
                                    --correction_methods {params.correction_methods} \
                                    --results_dir {params.results_dir} \
                                    --xenium_dir {params.xenium_dir} \
                                    --count_correction_dir {params.count_correction_dir} \
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


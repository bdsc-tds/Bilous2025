from pathlib import Path
import yaml
import itertools
import pandas as pd

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
xenium_cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
results_dir = Path(config['results_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
# probably only need to run for lognorm data
signal_integrity_thresholds = [0.5,0.7]
correction_methods = ['resolvi','resolvi_supervised'] + [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
normalisations = ['lognorm',]
layers = ['data',]
references = ['matched_reference_combo']
methods = ['rctd_class_aware']
levels = ['Level2.1']

n_neighbors = 10
n_permutations = 30
n_repeats = 5
top_n = 20
top_n_lr = 10
scoring = 'f1'
markers_mode = ['diffexpr']#,'common_markers'] #'/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/markers/cellmarker_cell_types_markers.json'

# resolvi params
num_samples = 30
mixture_k = 50

# needed to get unique cell types names for each level
# cell_types_palette = pd.read_csv(palette_dir / 'col_palette_cell_types_combo.csv')

out_files = []
for markers in markers_mode:
    for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
        if segmentation.stem == 'proseg_mode':
            continue
        for condition in (conditions := segmentation.iterdir()): 
            for panel in (panels := condition.iterdir()):
                for donor in (donors := panel.iterdir()):
                    for sample in (samples := donor.iterdir()):
                        for normalisation in normalisations:
                            for layer in layers:
                                for reference in references:
                                    for method in methods:
                                        for level in levels:
                                            for correction_method in correction_methods:

                                                k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                                                name = '/'.join(k)

                                                if 'proseg' in segmentation.stem:
                                                    k_proseg = ('proseg',condition.stem,panel.stem,donor.stem,sample.stem)
                                                    name_proseg = '/'.join(k_proseg)
                                                    sample_dir = xenium_dir / f'{name_proseg}/raw_results'
                                                else:
                                                    sample_dir = xenium_dir / f'{name}/normalised_results/outs'

                                                if correction_method == "resolvi":
                                                    name_corrected = f'{name}/{mixture_k=}/{num_samples=}/'
                                                elif correction_method == "resolvi_supervised":
                                                    name_corrected = f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/{mixture_k=}/{num_samples=}'
                                                elif "ovrlpy" in correction_method:
                                                    name_corrected = f'{name}'

                                                sample_corrected_counts_path = results_dir / f"{correction_method}/{name_corrected}/corrected_counts.h5"
                                                out_file_df_ctj_marker_genes = results_dir /  f'contamination_metrics_{markers}_corrected_counts/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_marker_genes.parquet'
                                                out_file_df_diffexpr = results_dir / f'contamination_metrics_{markers}_corrected_counts/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_diffexpr.parquet'
                                                out_file_df_markers_rank_significance_diffexpr = results_dir / f'contamination_metrics_{markers}_corrected_counts/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_markers_rank_significance_diffexpr.parquet'

                                                sample_normalised_counts = xenium_std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/{layer}.parquet'
                                                sample_idx = xenium_std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/cells.parquet'
                                                sample_annotation = xenium_cell_type_annotation_dir / f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet'

                                                if sample_corrected_counts_path.exists():

                                                    out_files.extend([
                                                        out_file_df_ctj_marker_genes,
                                                        out_file_df_diffexpr,
                                                        out_file_df_markers_rank_significance_diffexpr,
                                                        ])

                                                    rule:
                                                        name: f'contamination_metrics_{markers}_corrected_counts/{correction_method}/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                                                        input:
                                                            sample_corrected_counts_path=sample_corrected_counts_path,
                                                            sample_dir=sample_dir,
                                                            sample_normalised_counts=sample_normalised_counts,
                                                            sample_idx=sample_idx,
                                                            sample_annotation=sample_annotation,
                                                        output:
                                                            out_file_df_ctj_marker_genes=out_file_df_ctj_marker_genes,
                                                            out_file_df_diffexpr=out_file_df_diffexpr,
                                                            out_file_df_markers_rank_significance_diffexpr=out_file_df_markers_rank_significance_diffexpr,
                                                        params:
                                                            n_neighbors=n_neighbors,
                                                            n_permutations=n_permutations,
                                                            n_repeats=n_repeats,
                                                            top_n=top_n,
                                                            scoring=scoring,
                                                            markers=markers,
                                                        threads: 1
                                                        resources:
                                                            mem='30GB',
                                                            runtime='1h',
                                                        conda:
                                                            "spatial"
                                                        shell:
                                                            """
                                                            mkdir -p "$(dirname {output.out_file_df_diffexpr})"

                                                            python workflow/scripts/xenium/contamination_metrics_diffexpr_sample.py \
                                                                --sample_corrected_counts_path {input.sample_corrected_counts_path} \
                                                                --sample_dir {input.sample_dir} \
                                                                --sample_normalised_counts {input.sample_normalised_counts} \
                                                                --sample_idx {input.sample_idx} \
                                                                --sample_annotation {input.sample_annotation} \
                                                                --out_file_df_ctj_marker_genes {output.out_file_df_ctj_marker_genes} \
                                                                --out_file_df_diffexpr {output.out_file_df_diffexpr} \
                                                                --out_file_df_markers_rank_significance_diffexpr {output.out_file_df_markers_rank_significance_diffexpr} \
                                                                --n_neighbors {params.n_neighbors} \
                                                                --top_n {params.top_n} \
                                                                --scoring {params.scoring} \
                                                                --markers {params.markers} \

                                                            echo "DONE"
                                                            """


rule contamination_metrics_diffexpr_corrected_counts_all:
    input:
        out_files
    output:
        touch(results_dir / "contamination_metrics_{markers}_corrected_counts.done")



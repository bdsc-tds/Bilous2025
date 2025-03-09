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
normalisations = ['lognorm',]
layers = ['data',]
references = ['matched_reference_combo']
methods = ['rctd_class_aware']
levels = ['Level2']

n_neighbors = 10
n_permutations = 30
n_repeats = 5
top_n = 20
top_n_lr = 10
scoring = 'f1'
markers = 'diffexpr' #'/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/markers/cellmarker_cell_types_markers.json'

# needed to get unique cell types names for each level
# cell_types_palette = pd.read_csv(palette_dir / 'col_palette_cell_types_combo.csv')

out_files = []
for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    if segmentation.stem in ['proseg_mode','proseg_expected']:
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

                                        k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                                        name = '/'.join(k)

                                        if 'proseg' in segmentation.stem:
                                            k_proseg = ('proseg',condition.stem,panel.stem,donor.stem,sample.stem)
                                            name_proseg = '/'.join(k_proseg)
                                            sample_dir = xenium_dir / f'{name_proseg}/raw_results'
                                        else:
                                            sample_dir = xenium_dir / f'{name}/normalised_results/outs'

                                        sample_normalised_counts = xenium_std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/{layer}.parquet'
                                        sample_idx = xenium_std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/cells.parquet'
                                        sample_annotation = xenium_cell_type_annotation_dir / f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet'

                                        out_file_df_diffexpr = results_dir / f'contamination_metrics_diffexpr/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_diffexpr.parquet'
                                        out_file_df_markers_rank_significance_diffexpr = results_dir / f'contamination_metrics_diffexpr/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_markers_rank_significance_diffexpr.parquet'

                                        out_files.extend([
                                            out_file_df_diffexpr,
                                            out_file_df_markers_rank_significance_diffexpr,
                                            ])

                                        rule:
                                            name: f'contamination_metrics_diffexpr/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                                            input:
                                                sample_dir=sample_dir,
                                                sample_normalised_counts=sample_normalised_counts,
                                                sample_idx=sample_idx,
                                                sample_annotation=sample_annotation,
                                            output:
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
                                                    --sample_dir {input.sample_dir} \
                                                    --sample_normalised_counts {input.sample_normalised_counts} \
                                                    --sample_idx {input.sample_idx} \
                                                    --sample_annotation {input.sample_annotation} \
                                                    --out_file_df_diffexpr {output.out_file_df_diffexpr} \
                                                    --out_file_df_markers_rank_significance_diffexpr {output.out_file_df_markers_rank_significance_diffexpr} \
                                                    --n_neighbors {params.n_neighbors} \
                                                    --top_n {params.top_n} \
                                                    --scoring {params.scoring} \
                                                    --markers {params.markers} \

                                                echo "DONE"
                                                """


rule contamination_metrics_diffexpr_all:
    input:
        out_files
    output:
        touch(results_dir / "contamination_metrics_diffexpr.done")



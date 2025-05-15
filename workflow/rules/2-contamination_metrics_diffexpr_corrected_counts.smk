from pathlib import Path
import yaml
import itertools
import pandas as pd

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_count_correction_dir = Path(config['xenium_count_correction_dir'])
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
xenium_cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
results_dir = Path(config['results_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
# probably only need to run for lognorm data
signal_integrity_thresholds = [0.5,0.7]
correction_methods = ['split_fully_purified','resolvi','resolvi_supervised'] + [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
normalisations = ['lognorm',]
layers = ['data',]
references = ['matched_reference_combo']
methods = ['rctd_class_aware']
levels = ['Level2.1']


# params from pipeline config
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5


radius = 15
top_n = 20
markers_modes = ['diffexpr']#,'common_markers'] #'/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/markers/cellmarker_cell_types_markers.json'

genes_dict = {
    'all':[],# default: use all genes if empty list
    # 'Xenium_NSCLC_5k_lung_chromium_common_genes':pd.read_csv(config['markers_dir']+'Xenium_NSCLC_5k_lung_chromium_common_genes.csv')['gene'].tolist(),
    'Xenium_hLung_v1_metadata':pd.read_csv(config['markers_dir']+'Xenium_hLung_v1_metadata.csv')['Gene'].tolist(),
    # 'CHUV_IO_340_panel':pd.read_csv(config['markers_dir']+'CHUV_IO_340_panel.csv')['Gene ID'].tolist(),
    # 'Xenium_hBreast_v1_metadata':pd.read_csv(config['markers_dir']+'Xenium_hBreast_v1_metadata.csv')['Gene'].tolist()
}

# resolvi params
num_samples = 30
mixture_k = 50


out_files = []

for genes_name, genes in genes_dict.items():
    for markers_mode in markers_modes:
        for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
            if segmentation.stem == 'proseg_mode':
                continue
            for condition in (conditions := segmentation.iterdir()): 
                for panel in (panels := condition.iterdir()):
                    if genes_name != 'all' and panel.name != '5k':
                        continue
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
                                                    name_params = f"{markers_mode}_{radius=}_{top_n=}"

                                                    if 'proseg' in segmentation.stem:
                                                        k_proseg = ('proseg',condition.stem,panel.stem,donor.stem,sample.stem)
                                                        name_proseg = '/'.join(k_proseg)
                                                        sample_dir = xenium_dir / f'{name_proseg}/raw_results'
                                                    else:
                                                        sample_dir = xenium_dir / f'{name}/normalised_results/outs'

                                                    if correction_method == "split_fully_purified":
                                                        name_corrected = f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/split_fully_purified/'
                                                        sample_corrected_counts_path = xenium_count_correction_dir / f"{name_corrected}/corrected_counts.h5"

                                                    else:
                                                        if correction_method == "resolvi":
                                                            name_corrected = f'{name}/{mixture_k=}/{num_samples=}/'
                                                        elif correction_method == "resolvi_supervised":
                                                            name_corrected = f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/{mixture_k=}/{num_samples=}'
                                                        elif "ovrlpy" in correction_method:
                                                            name_corrected = f'{name}'

                                                        sample_corrected_counts_path = results_dir / f"{correction_method}/{name_corrected}/corrected_counts.h5"

                                                    output_dir = f'contamination_metrics_{name_params}/{genes_name}'

                                                    sample_normalised_counts = xenium_std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/{layer}.parquet'
                                                    sample_idx = xenium_std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/cells.parquet'
                                                    sample_annotation = xenium_cell_type_annotation_dir / f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet'
                                                    precomputed_ctj_markers = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_marker_genes.parquet'
                                                    precomputed_adata_obs = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_adata_obs.parquet'

                                                    out_file_df_ctj_marker_genes = results_dir /  f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_marker_genes.parquet'
                                                    out_file_df_diffexpr = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_diffexpr.parquet'
                                                    out_file_df_markers_rank_significance_diffexpr = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_markers_rank_significance_diffexpr.parquet'
                                                    out_file_summary_stats = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_summary_stats.json'
                                                    out_file_adata_obs = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_adata_obs.parquet'

                                                    if sample_corrected_counts_path.exists():

                                                        out_files.extend([
                                                            out_file_df_ctj_marker_genes,
                                                            out_file_df_diffexpr,
                                                            out_file_df_markers_rank_significance_diffexpr,
                                                            out_file_summary_stats,
                                                            out_file_adata_obs
                                                            ])

                                                        rule:
                                                            name: f'{output_dir}_corrected_counts/{correction_method}/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                                                            input:
                                                                sample_corrected_counts_path=sample_corrected_counts_path,
                                                                sample_dir=sample_dir,
                                                                sample_normalised_counts=sample_normalised_counts,
                                                                sample_idx=sample_idx,
                                                                sample_annotation=sample_annotation,
                                                                precomputed_ctj_markers=precomputed_ctj_markers,
                                                                precomputed_adata_obs=precomputed_adata_obs,
                                                            output:
                                                                out_file_df_ctj_marker_genes=out_file_df_ctj_marker_genes,
                                                                out_file_df_diffexpr=out_file_df_diffexpr,
                                                                out_file_df_markers_rank_significance_diffexpr=out_file_df_markers_rank_significance_diffexpr,
                                                                out_file_summary_stats=out_file_summary_stats,
                                                                out_file_adata_obs=out_file_adata_obs,
                                                            params:
                                                                radius=radius,
                                                                n_repeats=n_repeats,
                                                                top_n=top_n,
                                                                markers=markers_mode,
                                                                min_counts=min_counts,
                                                                min_features=min_features,
                                                                max_counts=max_counts,
                                                                max_features=max_features,
                                                                min_cells=min_cells,
                                                                genes=genes,
                                                            threads: 1
                                                            resources:
                                                                mem='50GB',
                                                                runtime='3h',
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
                                                                    --precomputed_ctj_markers {input.precomputed_ctj_markers} \
                                                                    --precomputed_adata_obs {input.precomputed_adata_obs} \
                                                                    --out_file_df_ctj_marker_genes {output.out_file_df_ctj_marker_genes} \
                                                                    --out_file_df_diffexpr {output.out_file_df_diffexpr} \
                                                                    --out_file_df_markers_rank_significance_diffexpr {output.out_file_df_markers_rank_significance_diffexpr} \
                                                                    --out_file_summary_stats {output.out_file_summary_stats} \
                                                                    --out_file_adata_obs {output.out_file_adata_obs} \
                                                                    --radius {params.radius} \
                                                                    --top_n {params.top_n} \
                                                                    --markers {params.markers} \
                                                                    --min_counts {params.min_counts} \
                                                                    --min_features {params.min_features} \
                                                                    --max_counts {params.max_counts} \
                                                                    --max_features {params.max_features} \
                                                                    --min_cells {params.min_cells} \
                                                                    --genes {params.genes} \

                                                                echo "DONE"
                                                                """


rule contamination_metrics_diffexpr_corrected_counts_all:
    input:
        out_files
    output:
        touch(results_dir / f"{output_dir}_corrected_counts.done")



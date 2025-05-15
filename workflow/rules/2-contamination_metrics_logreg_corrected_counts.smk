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
n_permutations = 30
n_splits= 5
top_n = 20
scoring = 'precision'
cv_mode = 'spatial'
markers_modes = ['diffexpr']#,'common_markers'] #'/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/markers/cellmarker_cell_types_markers.json'
max_n_cells = 50_000

# resolvi params
num_samples = 30
mixture_k = 50

genes_dict = {
    'all':[],# default: use all genes if empty list
    # 'Xenium_NSCLC_5k_lung_chromium_common_genes':pd.read_csv(config['markers_dir']+'Xenium_NSCLC_5k_lung_chromium_common_genes.csv')['gene'].tolist(),
    'Xenium_hLung_v1_metadata':pd.read_csv(config['markers_dir']+'Xenium_hLung_v1_metadata.csv')['Gene'].tolist(),
    # 'CHUV_IO_340_panel':pd.read_csv(config['markers_dir']+'CHUV_IO_340_panel.csv')['Gene ID'].tolist(),
    # 'Xenium_hBreast_v1_metadata':pd.read_csv(config['markers_dir']+'Xenium_hBreast_v1_metadata.csv')['Gene'].tolist()
}


out_files = []
for train_mode in train_modes:
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
                                                        name_params_diffexpr = f"{markers_mode}_{radius=}_{top_n=}"
                                                        name_params = f"{markers_mode}_{radius=}_{n_permutations=}_{n_splits=}_{top_n=}_{scoring}_{cv_mode}_{train_mode}"

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

                                                        output_dir_diffexpr = f'contamination_metrics_{name_params_diffexpr}/{genes_name}'
                                                        output_dir = f'contamination_metrics_{name_params}_logreg/{genes_name}'

                                                        sample_normalised_counts = xenium_std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/{layer}.parquet'
                                                        sample_idx = xenium_std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/cells.parquet'
                                                        sample_annotation = xenium_cell_type_annotation_dir / f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet'
                                                        precomputed_ctj_markers = results_dir / f'{output_dir_diffexpr}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_marker_genes.parquet'
                                                        precomputed_adata_obs = results_dir / f'{output_dir_diffexpr}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_adata_obs.parquet'
                                                        
                                                        out_file_df_permutations_logreg = results_dir / f'{output_dir}/{correction_method}/{name_corrected}/{normalisation}/{layer}_{reference}_{method}_{level}_permutations_logreg.parquet'
                                                        out_file_df_importances_logreg = results_dir / f'{output_dir}/{correction_method}/{name_corrected}//{normalisation}/{layer}_{reference}_{method}_{level}_importances_logreg.parquet'
                                                        out_file_df_markers_rank_significance_logreg = results_dir / f'{output_dir}/{correction_method}/{name_corrected}//{normalisation}/{layer}_{reference}_{method}_{level}_markers_rank_significance_logreg.parquet'

                                                        if sample_corrected_counts_path.exists():

                                                            out_files.extend([
                                                                out_file_df_permutations_logreg,
                                                                out_file_df_importances_logreg,
                                                                out_file_df_markers_rank_significance_logreg,
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
                                                                    out_file_df_permutations_logreg=out_file_df_permutations_logreg,
                                                                    out_file_df_importances_logreg=out_file_df_importances_logreg,
                                                                    out_file_df_markers_rank_significance_logreg=out_file_df_markers_rank_significance_logreg,
                                                                params:
                                                                    radius=radius,
                                                                    n_splits=n_splits,
                                                                    n_permutations=n_permutations,
                                                                    n_repeats=n_repeats,
                                                                    cv_mode=cv_mode,
                                                                    top_n=top_n,
                                                                    scoring=scoring,
                                                                    markers=markers_mode,
                                                                    max_n_cells=max_n_cells,
                                                                    min_counts=min_counts,
                                                                    min_features=min_features,
                                                                    max_counts=max_counts,
                                                                    max_features=max_features,
                                                                    min_cells=min_cells,
                                                                    genes=genes,
                                                                    train_mode=train_mode,
                                                                threads: 1
                                                                resources:
                                                                    mem='50GB',
                                                                    runtime='2d',
                                                                conda:
                                                                    "spatial"
                                                                shell:
                                                                    """
                                                                    mkdir -p "$(dirname {output.out_file_df_permutations_logreg})"

                                                                    python workflow/scripts/xenium/contamination_metrics_logreg_sample.py \
                                                                        --sample_corrected_counts_path {input.sample_corrected_counts_path} \
                                                                        --sample_dir {input.sample_dir} \
                                                                        --sample_normalised_counts {input.sample_normalised_counts} \
                                                                        --sample_idx {input.sample_idx} \
                                                                        --sample_annotation {input.sample_annotation} \
                                                                        --precomputed_ctj_markers {input.precomputed_ctj_markers} \
                                                                        --precomputed_adata_obs {input.precomputed_adata_obs} \
                                                                        --out_file_df_permutations_logreg {output.out_file_df_permutations_logreg} \
                                                                        --out_file_df_importances_logreg {output.out_file_df_importances_logreg} \
                                                                        --out_file_df_markers_rank_significance_logreg {output.out_file_df_markers_rank_significance_logreg} \
                                                                        --radius {params.radius} \
                                                                        --n_splits {params.n_splits} \
                                                                        --n_permutations {params.n_permutations} \
                                                                        --n_repeats {params.n_repeats} \
                                                                        --cv_mode {params.cv_mode} \
                                                                        --top_n {params.top_n} \
                                                                        --scoring {params.scoring} \
                                                                        --markers {params.markers} \
                                                                        --max_n_cells {params.max_n_cells} \
                                                                        --min_counts {params.min_counts} \
                                                                        --min_features {params.min_features} \
                                                                        --max_counts {params.max_counts} \
                                                                        --max_features {params.max_features} \
                                                                        --min_cells {params.min_cells} \
                                                                        --genes {params.genes} \
                                                                        --train_mode {params.train_mode} \

                                                                    echo "DONE"
                                                                    """


rule contamination_metrics_logreg_corrected_counts_all:
    input:
        out_files
    output:
        touch(results_dir / f"{output_dir}_corrected_counts.done")



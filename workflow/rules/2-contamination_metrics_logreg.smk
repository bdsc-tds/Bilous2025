train_modes = ['multivariate','univariate']
params_product = list(product(normalisations, layers, references, methods, levels))

out_files = []
for train_mode in train_modes:
    for genes_name, genes in genes_dict.items():
        for markers_mode in markers_modes:
            for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
                if segmentation.stem in ['proseg_mode','bats_normalised','bats_expected']:
                    continue
                for condition in (conditions := segmentation.iterdir()): 
                    for panel in (panels := condition.iterdir()):
                        if genes_name != 'all' and panel.name != '5k':
                            continue
                        for donor in (donors := panel.iterdir()):
                            for sample in (samples := donor.iterdir()):
                                for normalisation, layer, reference, method, level in params_product:
                                    if level not in CONDITIONS_LEVELS[condition.stem]:
                                        continue
                                    if reference not in CONDITIONS_REFERENCES[condition.stem]:
                                        continue

                                    k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem)
                                    name = '/'.join(k)

                                    name_params_diffexpr = f"{markers_mode}_{radius=}_{top_n=}"
                                    name_params = f"{markers_mode}_{radius=}_{n_permutations=}_{n_splits=}_{top_n=}_{scoring}_{cv_mode}_{train_mode}"

                                    if 'proseg' in segmentation.stem:
                                        k_proseg = ('proseg',condition.stem,panel.stem,donor.stem,sample.stem)
                                        name_proseg = '/'.join(k_proseg)
                                        sample_dir = xenium_processed_data_dir / f'{name_proseg}/raw_results'
                                    else:
                                        sample_dir = xenium_processed_data_dir / f'{name}/normalised_results/outs'

                                    output_dir_diffexpr = f'contamination_metrics_{name_params_diffexpr}/{genes_name}'
                                    output_dir = f'contamination_metrics_{name_params}_logreg/{genes_name}'

                                    sample_normalised_counts = std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/{layer}.parquet'
                                    sample_idx = std_seurat_analysis_dir / f'{name}/{normalisation}/normalised_counts/cells.parquet'
                                    sample_annotation = cell_type_annotation_dir / f'{name}/{normalisation}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet'
                                    precomputed_ctj_markers = results_dir / f'{output_dir_diffexpr}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_marker_genes.parquet'
                                    precomputed_adata_obs = results_dir / f'{output_dir_diffexpr}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_adata_obs.parquet'
                                    
                                    out_file_df_permutations_logreg = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_permutations_logreg.parquet'
                                    out_file_df_importances_logreg = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_importances_logreg.parquet'
                                    out_file_df_markers_rank_significance_logreg = results_dir / f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}_markers_rank_significance_logreg.parquet'

                                    out_files.extend([
                                        out_file_df_permutations_logreg,
                                        out_file_df_importances_logreg,
                                        out_file_df_markers_rank_significance_logreg,
                                        ])

                                    rule:
                                        name: f'{output_dir}/raw/{name}/{normalisation}/{layer}_{reference}_{method}_{level}'
                                        input:
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
                                                --genes {params.genes} \
                                                --train_mode {params.train_mode} \

                                            echo "DONE"
                                            """


rule contamination_metrics_logreg_all:
    input:
        out_files
    output:
        touch(results_dir / f"{output_dir}.done")



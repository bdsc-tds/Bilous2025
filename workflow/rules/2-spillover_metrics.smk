from pathlib import Path
import yaml

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
xenium_cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
results_dir = Path(config['results_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
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

# cell_type_palette = pd.read_csv(palette_dir / 'col_palette_cell_types_combo.csv')

out_files = []
for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    for normalisation in normalisations:
                        for layer in layers:
                            for reference in references:
                                for method in methods:
                                    for level in levels:
                                        k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,normalisation)
                                        name = '/'.join(k)

                                        sample_normalised_counts = xenium_std_seurat_analysis_dir / f'{name}/normalised_counts/{layer}.parquet'
                                        sample_idx = xenium_std_seurat_analysis_dir / f'{name}/normalised_counts/cells.parquet'
                                        sample_annotation = xenium_cell_type_annotation_dir / f'{name}/reference_based/{reference}/{method}/{level}/single_cell/labels.parquet'

                                        if not sample_annotation_dir.exists():
                                            continue
                                            
                                        out_file_df_permutations_logreg = results_dir / f'contamination_metrics/{name}/{layer}_permutations_logreg.parquet'
                                        out_file_df_importances_logreg = results_dir / f'contamination_metrics/{name}/{layer}_importances_logreg.parquet'
                                        out_file_df_diffexpr = results_dir / f'contamination_metrics/{name}/{layer}_diffexpr.parquet'
                                        out_file_df_markers_rank_significance_logreg = results_dir / f'contamination_metrics/{name}/{layer}_markers_rank_significance_logreg.json'
                                        out_file_df_markers_rank_significance_diffexpr = results_dir / f'contamination_metrics/{name}/{layer}_markers_rank_significance_diffexpr.json'
                                        out_dir_liana_lrdata = results_dir / f'contamination_metrics/{name}/{layer}_liana_lrdata_folder'

                                        out_files.extend([
                                            out_file_df_permutations_logreg,
                                            out_file_df_importances_logreg,
                                            out_file_df_diffexpr,
                                            out_file_df_markers_rank_significance_logreg,
                                            out_file_df_markers_rank_significance_diffexpr,
                                            out_dir_liana_lrdata])

                                        rule:
                                            name: f'contamination_metrics/{name}/{layer}'
                                            input:
                                                sample=sample,
                                                sample_normalised_counts=sample_normalised_counts,
                                                sample_idx=sample_idx,
                                                sample_annotation=sample_annotation,
                                            output:
                                                out_file_df_permutations_logreg=out_file_df_permutations_logreg,
                                                out_file_df_importances_logreg=out_file_df_importances_logreg,
                                                out_file_df_diffexpr=out_file_df_diffexpr,
                                                out_file_df_markers_rank_significance_logreg=out_file_df_markers_rank_significance_logreg,
                                                out_file_df_markers_rank_significance_diffexpr=out_file_df_markers_rank_significance_diffexpr,
                                                out_dir_liana_lrdata=out_dir_liana_lrdata,
                                            params:
                                                n_neighbors=n_neighbors,
                                                n_permutations=n_permutations,
                                                n_repeats=n_repeats,
                                                top_n=top_n,
                                                top_n_lr=top_n_lr,
                                                cti=cti,
                                                ctj=ctj,
                                                scoring=scoring,
                                                markers=markers,
                                            threads: 1
                                            resources:
                                                mem='30GB',
                                                runtime='15m',
                                            conda:
                                                "spatial"
                                            shell:
                                                """
                                                mkdir -p "$(dirname {output.out_file_df_permutations_logreg})"

                                                python workflow/scripts/xenium/contamination_metrics_sample.py \
                                                    --sample {input.sample} \
                                                    --sample_normalised_counts {input.sample_normalised_counts} \
                                                    --sample_idx {input.sample_idx} \
                                                    --sample_annotation {input.sample_annotation} \
                                                    --out_file_df_permutations_logreg {output.out_file_df_permutations_logreg} \
                                                    --out_file_df_importances_logreg {output.out_file_df_importances_logreg} \
                                                    --out_file_df_diffexpr {output.out_file_df_diffexpr} \
                                                    --out_file_df_markers_rank_significance_logreg {output.out_file_df_markers_rank_significance_logreg} \
                                                    --out_file_df_markers_rank_significance_diffexpr {output.out_file_df_markers_rank_significance_diffexpr} \
                                                    --out_dir_liana_lrdata {output.out_dir_liana_lrdata} \
                                                    --n_neighbors {params.n_neighbors} \
                                                    --n_permutations {params.n_permutations} \
                                                    --n_repeats {params.n_repeats} \
                                                    --top_n {params.top_n} \
                                                    --top_n_lr {params.top_n_lr} \
                                                    --cti {params.cti} \
                                                    --ctj {params.ctj} \
                                                    --scoring {params.scoring} \
                                                    --markers {params.markers} \

                                                echo "DONE"
                                                """


rule spillover_metrics_all:
    input:
        out_files
    output:
        touch(results_dir / "spillover_metrics.done")



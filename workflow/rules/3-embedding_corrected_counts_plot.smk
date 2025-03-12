# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
xenium_cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
n_comps = config['umap_n_comps']
n_neighbors = config['umap_n_neighbors']
min_dist = config['umap_min_dist']
metric = config['umap_metric']
s=0.5
alpha=0.5

cell_type_palette = palette_dir / 'col_palette_cell_types_combo.csv'
panel_palette = palette_dir / 'col_palette_panel.csv'
sample_palette = palette_dir / 'col_palette_sample.csv'

signal_integrity_thresholds = [0.5,0.7]
correction_methods = ['resolvi','resolvi_supervised'] + [f'ovrlpy_correction_{signal_integrity_threshold=}' for signal_integrity_threshold in signal_integrity_thresholds]
normalisations = ['lognorm','sctransform']
layers = ['data','scale_data']
references = ['matched_reference_combo','external_reference']
methods = ['rctd_class_aware']
colors = ['sample','Level2.1']#['Level1','Level2','Level3','Level4','panel','sample',] # condition and sample as color to plot added here in addition to levels
extension = 'png'

out_files_panel = []

for correction_method in correction_methods:
    for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
        if segmentation.stem == 'proseg_mode':
            continue
        for condition in (conditions := segmentation.iterdir()): 
            for panel in (panels := condition.iterdir()):
                # for normalisation in normalisations:
                #     for layer in layers: 
                        for reference in references:
                            for method in methods:
                                for color in colors:

                                    if color == 'Level2.1' and reference == 'external_reference':
                                        continue

                                    # input embedding file (doesn't depend on ref,method or color loops but more readable to have here)
                                    k = (segmentation.stem,condition.stem,panel.stem)#,normalisation)
                                    name = '/'.join(k)                                                      # {layer}_
                                    embed_file = results_dir / f'{correction_method}_embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet'

                                    # no need to plot panel for panel level UMAPs
                                    if color == 'panel':
                                        continue
                                    
                                    # no need to plot sample coloring for every param combination
                                    if color == 'sample' and (reference != references[0] or method != methods[0]):
                                        continue
                                                                                                                                                                    # _{layer}
                                    out_file = figures_dir / f"{correction_method}_embed_panel/{name}/umap_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}_{reference}_{method}_{color}.{extension}"
                                    out_files_panel.append(out_file)

                                    rule:
                                        name: f'{correction_method}_embed_panel_plot/{name}/umap_{reference}_{method}_{color}'#_{layer}'
                                        input:
                                            panel=panel,
                                            embed_file=embed_file,
                                        output:
                                            out_file=out_file,
                                        params:
                                            cell_type_annotation_dir=xenium_cell_type_annotation_dir,
                                            normalisation=normalisation,
                                            reference=reference,
                                            method=method,
                                            color=color,
                                            cell_type_palette=cell_type_palette,
                                            panel_palette=panel_palette,
                                            sample_palette=sample_palette,
                                            s=s,
                                            alpha=alpha,
                                        threads: 1
                                        resources:
                                            mem='30GB',
                                            runtime='10m',
                                        conda:
                                            "spatial"
                                        shell:
                                            """
                                            mkdir -p "$(dirname {output.out_file})"

                                            python workflow/scripts/xenium/embed_panel_plot.py \
                                            --panel {input.panel} \
                                            --embed_file {input.embed_file} \
                                            --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                                            --normalisation {params.normalisation} \
                                            --reference {params.reference} \
                                            --method {params.method} \
                                            --color {params.color} \
                                            --out_file {output.out_file} \
                                            --cell_type_palette {params.cell_type_palette} \
                                            --panel_palette {params.panel_palette} \
                                            --sample_palette {params.sample_palette} \
                                            --s {params.s} \
                                            --alpha {params.alpha} \
                                            
                                            echo "DONE"
                                            """


rule embed_panel_corrected_counts_plot_all:
    input:
        out_files_panel

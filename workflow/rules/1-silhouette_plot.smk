from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
silhouette_dir = Path(config['results_dir']) / 'silhouette'
figures_dir = Path(config['figures_dir'])
palette_dir = Path(config['xenium_metadata_dir'])

# Params
metric='euclidean'
segmentation = list(xenium_dir.iterdir())[0] # take arbitrary segmentation just to loop over folders
segmentation_palette = palette_dir / 'col_palette_segmentation.csv'
normalisations = ['lognorm','sctransform']
layers = ['data','scale_data']
references = ['matched_reference_combo','external_reference']
methods = ['rctd_class_aware']
levels = ['Level2'] 
extension = 'png'
dpi = 100

out_files_panel = []
for condition in (conditions := segmentation.iterdir()): 
    for panel in (panels := condition.iterdir()):
        for normalisation in normalisations:
            for layer in layers:
                for reference in references:
                    for method in methods:
                        for level in levels:

                            k = (condition.stem,panel.stem,normalisation)
                            name = '/'.join(k)

                            out_file = figures_dir / f'silhouette_plot_panel/{name}/silhouette_{layer}_{metric}_{reference}_{method}_{level}.{extension}'
                            out_files_panel.append(out_file)

                            rule:
                                name: f'silhouette_plot_panel/{name}/silhouette_{layer}_{metric}_{reference}_{method}_{level}'
                                input:
                                    silhouette_is_done=results_dir / "silhouette.done",                            
                                output:
                                    out_file=out_file,
                                params:
                                    silhouette_dir=silhouette_dir,
                                    plot_condition=condition.stem,
                                    plot_panel=panel.stem,
                                    segmentation_palette=segmentation_palette,
                                    normalisation=normalisation,
                                    layer=layer,
                                    reference=reference,
                                    method=method,
                                    level=level,
                                    dpi=dpi,
                                threads: 1
                                resources:
                                    mem='5GB',
                                    runtime='5m',
                                conda:
                                    "spatial"
                                shell:
                                    """
                                    mkdir -p "$(dirname {output.out_file})"

                                    python workflow/scripts/xenium/silhouette_panel_plot.py \
                                    --silhouette_dir {params.silhouette_dir} \
                                    --plot_condition {params.plot_condition} \
                                    --plot_panel {params.plot_panel} \
                                    --out_file {output.out_file} \
                                    --segmentation_palette {params.segmentation_palette} \
                                    --normalisation {params.normalisation} \
                                    --layer {params.layer} \
                                    --reference {params.reference} \
                                    --method {params.method} \
                                    --level {params.level} \
                                    --dpi {params.dpi}

                                    echo "DONE"
                                    """





out_files_condition = []

for condition in (conditions := list(segmentation.iterdir())+[Path('all')]): 
        for normalisation in normalisations:
            for layer in layers:
                for reference in references:
                    for method in methods:
                        for level in levels:

                            k = (condition.stem,normalisation)
                            name = '/'.join(k)

                            out_file = figures_dir / f'silhouette_plot_condition/{name}/silhouette_{layer}_{metric}_{reference}_{method}_{level}.{extension}'
                            out_files_condition.append(out_file)

                            rule:
                                name: f'silhouette_plot_condition/{name}/silhouette_{layer}_{metric}_{reference}_{method}_{level}'
                                input:
                                    silhouette_is_done=results_dir / "silhouette.done",                            
                                output:
                                    out_file=out_file,
                                params:
                                    silhouette_dir=silhouette_dir,
                                    plot_condition=condition.stem,
                                    segmentation_palette=segmentation_palette,
                                    normalisation=normalisation,
                                    layer=layer,
                                    reference=reference,
                                    method=method,
                                    level=level,
                                    dpi=dpi,
                                threads: 1
                                resources:
                                    mem='5GB',
                                    runtime='5m',
                                conda:
                                    "spatial"
                                shell:
                                    """
                                    mkdir -p "$(dirname {output.out_file})"

                                    python workflow/scripts/xenium/silhouette_condition_plot.py \
                                    --silhouette_dir {params.silhouette_dir} \
                                    --plot_condition {params.plot_condition} \
                                    --out_file {output.out_file} \
                                    --segmentation_palette {params.segmentation_palette} \
                                    --normalisation {params.normalisation} \
                                    --layer {params.layer} \
                                    --reference {params.reference} \
                                    --method {params.method} \
                                    --level {params.level} \
                                    --dpi {params.dpi}

                                    echo "DONE"
                                    """


rule silhouette_plot_panel_all:
    input:
        out_files_panel


rule silhouette_plot_condition_all:
    input:
        out_files_condition
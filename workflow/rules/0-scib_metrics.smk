# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])

# Params
n_comps = 50
max_n_cells = 1_000_000

normalisations = ['lognorm','sctransform']
layers = ['data','scale_data']
references = ['matched_reference_combo']#,'external_reference']
methods = ['rctd_class_aware']
levels = ['Level2'] # condition and sample as color to plot added here in addition to levels

out_files_panel = []

for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for normalisation in normalisations:
                for layer in layers: 
                    for reference in references:
                        for method in methods:
                            for level in levels:

                                # input embedding file (doesn't depend on ref,method or color loops but more readable to have here)
                                k = (segmentation.stem,condition.stem,panel.stem,normalisation)
                                name = '/'.join(k)

                                out_file = results_dir / f"scib_metrics_panel/{name}/scib_metrics_{layer}_{reference}_{method}_{level}_{n_comps=}_{max_n_cells=}.parquet"
                                out_files_panel.append(out_file)

                                rule:
                                    name: f'scib_metrics_panel/{name}/scib_metrics_{layer}_{reference}_{method}_{level}'
                                    input:
                                        panel=panel,
                                    output:
                                        out_file=out_file,
                                    params:
                                        cell_type_annotation_dir=cell_type_annotation_dir,
                                        normalisation=normalisation,
                                        layer=layer,
                                        reference=reference,
                                        method=method,
                                        level=level,
                                        n_comps=n_comps,
                                        max_n_cells=max_n_cells,
                                    threads: 1
                                    resources:
                                        mem='80GB',
                                        runtime='5h',
                                    conda:
                                        "spatial"
                                    shell:
                                        """
                                        mkdir -p "$(dirname {output.out_file})"

                                        python workflow/scripts/xenium/scib_metrics_panel.py \
                                        --panel {input.panel} \
                                        --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                                        --out_file {output.out_file} \
                                        --normalisation {params.normalisation} \
                                        --layer {params.layer} \
                                        --reference {params.reference} \
                                        --method {params.method} \
                                        --level {params.level} \
                                        --n_comps {params.n_comps} \
                                        --max_n_cells {params.max_n_cells} \
                                        
                                        echo "DONE"
                                        """



rule scib_metrics_panel_all:
    input:
        out_files_panel

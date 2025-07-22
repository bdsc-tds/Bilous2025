out_files_panel = []
for segmentation in (segmentations := xenium_processed_data_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for normalisation in normalisations:
                for layer in layers:
                    for reference in references:
                        for method in methods:
                            for level in levels:

                                k = (segmentation.stem,condition.stem,panel.stem,normalisation)
                                name = '/'.join(k)

                                panel_silhouette = results_dir / f'silhouette/{name}'
                                out_file = figures_dir / f'silhouette_plot_panel/{name}/silhouette_{layer}_{metric}_{reference}_{method}_{level}.{extension}'
                                out_files_panel.append(out_file)

                                rule:
                                    name: f'silhouette_plot_panel/{name}/silhouette_{layer}_{metric}_{reference}_{method}_{level}'
                                    input:
                                        silhouette_is_done=results_dir / "silhouette.done",                            
                                    output:
                                        out_file=out_file,
                                    params:
                                        panel=panel_silhouette,
                                        segmentation_palette=segmentation_palette,
                                        normalisation=normalisation,
                                        layer=layer,
                                        reference=reference,
                                        method=method,
                                        level=level,
                                        metric=metric,
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
                                        --panel {params.panel} \
                                        --out_file {output.out_file} \
                                        --segmentation_palette {params.segmentation_palette} \
                                        --normalisation {params.normalisation} \
                                        --layer {params.layer} \
                                        --reference {params.reference} \
                                        --method {params.method} \
                                        --level {params.level} \
                                        --metric {params.metric} \

                                        echo "DONE"
                                        """





out_files_condition = []
for segmentation in (segmentations := xenium_processed_data_dir.iterdir()):
    if segmentation.stem == 'proseg_v1':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for normalisation in normalisations:
            for reference in references:
                for method in methods:
                    for level in levels:

                        k = (segmentation.stem,condition.stem,normalisation)
                        name = '/'.join(k)
                        condition_silhouette = results_dir / f'silhouette/{name}'

                        out_file = figures_dir / f'silhouette_plot_condition/{name}/silhouette_{layer}_{metric}_{reference}_{method}_{level}.{extension}'
                        out_files_condition.append(out_file)

                        rule:
                            name: f'silhouette_plot_condition/{name}/silhouette_{layer}_{metric}_{reference}_{method}_{level}'
                            input:
                                condition=condition_silhouette,
                            output:
                                out_file=out_file,
                            params:
                                segmentation_palette=segmentation_palette,
                                normalisation=normalisation,
                                layer=layer,
                                reference=reference,
                                method=method,
                                level=level,
                                metric=metric,
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
                                --condition {input.condition} \
                                --out_file {output.out_file} \
                                --segmentation_palette {params.segmentation_palette} \
                                --normalisation {params.normalisation} \
                                --layer {params.layer} \
                                --reference {params.reference} \
                                --method {params.method} \
                                --level {params.level} \
                                --metric {params.metric} \

                                echo "DONE"
                                """


rule silhouette_plot_panel_all:
    input:
        out_files_panel


rule silhouette_plot_condition_all:
    input:
        out_files_condition
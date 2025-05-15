from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])

# params from pipeline config
radius = 50
b_cell_label = 'B cell'
t_cell_label = 'T cell'
min_perc_b = 0.0 # disable min % filter
min_perc_t = 0.0 # disable min % filter
fold_change_threshold = 1.5
min_b = 5
min_t = 5
min_tls_size = 20

# params supervised
normalisation = 'lognorm'
mode = 'reference_based'
references = ['matched_reference_combo']
methods = ['rctd_class_aware']
levels = ['Level2.1']

out_files = []
for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    if segmentation.stem == 'proseg_mode':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):

                    for reference in references:
                        for method in methods:
                            for level in levels:
                                
                                if segmentation.stem == 'proseg_expected':
                                    name_sample =  '/'.join(('proseg',condition.stem,panel.stem,donor.stem,sample.stem,))
                                    path = xenium_dir / f'{name_sample}/raw_results'
                                else:
                                    name_sample =  '/'.join((segmentation.stem.replace('proseg_mode','proseg'),condition.stem,panel.stem,donor.stem,sample.stem,))
                                    path = xenium_dir / f'{name_sample}/normalised_results/outs'

                                k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,
                                     normalisation,mode,reference,method,level)
                                name = '/'.join(k)
                                cell_type_labels = cell_type_annotation_dir / name / f"single_cell/labels.parquet"

                                if path.exists():

                                    out_file =  results_dir / f'detect_tls/{name}/obs_tls.parquet'
                                    out_files.append(out_file)

                                    rule:
                                        name: f'detect_tls/{name}'
                                        input:
                                            path=path,
                                        output:
                                            out_file=out_file,
                                        params:
                                            radius=radius,
                                            b_cell_label=b_cell_label,
                                            t_cell_label=t_cell_label,
                                            min_perc_b=min_perc_b,
                                            min_perc_t=min_perc_t,
                                            fold_change_threshold=fold_change_threshold,
                                            min_b=min_b,
                                            min_t=min_t,
                                            min_tls_size=min_tls_size,
                                            cell_type_labels=cell_type_labels,
                                        threads: 1
                                        resources:
                                            mem='50GB',# if panel.stem == '5k' else '10GB',
                                            runtime='15m',
                                        conda:
                                            "spatial"
                                        shell:
                                            """
                                            mkdir -p "$(dirname {output.out_file})"

                                            python workflow/scripts/xenium/detect_tls_sample.py \
                                            --path {input.path} \
                                            --out_file {output.out_file} \
                                            --radius {params.radius} \
                                            --b_cell_label "{params.b_cell_label}" \
                                            --t_cell_label "{params.t_cell_label}" \
                                            --min_perc_b {params.min_perc_b} \
                                            --min_perc_t {params.min_perc_t} \
                                            --fold_change_threshold {params.fold_change_threshold} \
                                            --min_b {params.min_b} \
                                            --min_t {params.min_t} \
                                            --min_tls_size {params.min_tls_size} \
                                            --cell_type_labels {params.cell_type_labels} 
                                            
                                            echo "DONE"
                                            """



rule detect_tls_all:
    input:
        out_files
    output:
        touch(results_dir / "detect_tls.done")
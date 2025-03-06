from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
xenium_std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])

# params from pipeline config
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# params
num_samples = 30
batch_size = 1000
mixture_k = 50

out_files_training = []
for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    
                    name_sample =  '/'.join((segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,))
                    if segmentation.stem == 'proseg_expected':
                        path = xenium_dir / f'{name_sample}/raw_results'
                    else:
                        path = xenium_dir / f'{name_sample}/normalised_results/outs'

                    k = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,f'{mixture_k=}')
                    name = '/'.join(k)

                    if path.exists():

                        out_dir_resolvi_model = results_dir / f'resolvi/{name}/model/'
                        out_file_resolvi_model = out_dir_resolvi_model / 'model.pt'
                        out_files_training.append(out_file_resolvi_model)

                        rule:
                            name: f'resolvi/{name}'
                            input:
                                path=path,
                            output:
                                out_file_resolvi_model = out_file_resolvi_model
                            params:
                                out_dir_resolvi_model=out_dir_resolvi_model,
                                min_counts=min_counts,
                                min_features=min_features,
                                max_counts=max_counts,
                                max_features=max_features,
                                min_cells=min_cells,
                                mixture_k=mixture_k,
                            threads: 1
                            resources:
                                mem='80GB',# if panel.stem == '5k' else '10GB',
                                runtime='3h',
                                slurm_partition = "gpu",
                                slurm_extra = '--gres=gpu:1',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {params.out_dir_resolvi_model})"

                                python workflow/scripts/xenium/resolvi_sample_training.py \
                                --path {input.path} \
                                --out_dir_resolvi_model {params.out_dir_resolvi_model} \
                                --min_counts {params.min_counts} \
                                --min_features {params.min_features} \
                                --max_counts {params.max_counts} \
                                --max_features {params.max_features} \
                                --min_cells {params.min_cells} \
                                --mixture_k {params.mixture_k} \
                                
                                echo "DONE"
                                """


out_files_inference = []
for segmentation in (segmentations := xenium_std_seurat_analysis_dir.iterdir()):
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for donor in (donors := panel.iterdir()):
                for sample in (samples := donor.iterdir()):
                    
                    k_sample = (segmentation.stem,condition.stem,panel.stem,donor.stem,sample.stem,)
                    name_sample =  '/'.join(k_sample)
                    if segmentation.stem == 'proseg_expected':
                        path = xenium_dir / f'{name_sample}/raw_results'
                    else:
                        path = xenium_dir / f'{name_sample}/normalised_results/outs'
                    
                    k_model = k_sample+(f'{mixture_k=}',)
                    k = k_model + (f'{num_samples=}',)
                    name_model = '/'.join(k_model)
                    name = '/'.join(k)

                    if path.exists():

                        dir_resolvi_model = results_dir / f'resolvi/{name_model}/model/'
                        out_file_resolvi_corrected_counts = results_dir / f'resolvi/{name}/corrected_counts.h5'
                        out_file_resolvi_proportions = results_dir/ f'resolvi/{name}/proportions.parquet'
                        out_files_inference.extend([out_file_resolvi_corrected_counts,out_file_resolvi_proportions])

                        rule:
                            name: f'resolvi_inference/{name}'
                            input:
                                path=path,
                                file_resolvi_model = dir_resolvi_model / 'model.pt'
                            output:
                                out_file_resolvi_corrected_counts=out_file_resolvi_corrected_counts,
                                out_file_resolvi_proportions=out_file_resolvi_proportions,
                            params:
                                dir_resolvi_model=dir_resolvi_model,
                                min_counts=min_counts,
                                min_features=min_features,
                                max_counts=max_counts,
                                max_features=max_features,
                                min_cells=min_cells,
                                num_samples=num_samples,
                                batch_size=batch_size,
                                mixture_k=mixture_k,
                            threads: 1
                            resources:
                                mem='200GB',# if panel.stem == '5k' else '10GB',
                                runtime='6h',
                                slurm_partition = "gpu",
                                slurm_extra = '--gres=gpu:1',
                            conda:
                                "spatial"
                            shell:
                                """
                                mkdir -p "$(dirname {output.out_file_resolvi_corrected_counts})"

                                python workflow/scripts/xenium/resolvi_sample_inference.py \
                                --path {input.path} \
                                --dir_resolvi_model {params.dir_resolvi_model} \
                                --out_file_resolvi_corrected_counts {output.out_file_resolvi_corrected_counts} \
                                --out_file_resolvi_proportions {output.out_file_resolvi_proportions} \
                                --min_counts {params.min_counts} \
                                --min_features {params.min_features} \
                                --max_counts {params.max_counts} \
                                --max_features {params.max_features} \
                                --min_cells {params.min_cells} \
                                --num_samples {params.num_samples} \
                                --batch_size {params.batch_size} \
                                
                                echo "DONE"
                                """


rule resolvi_training_all:
    input:
        out_files_training

rule resolvi_inference_all:
    input:
        out_files_inference
    output:
        touch(results_dir / "resolvi.done")
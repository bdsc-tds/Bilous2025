from pathlib import Path

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir']) 
scrnaseq_processed_data_dir = Path(config['scrnaseq_processed_data_dir'])
seurat_to_h5_dir = results_dir / 'seurat_to_h5'

# Params
genes = 'conditions'
methods = ['conditional','jaccard']#,'pearson','spearman']
target_counts = [30,50,200]
out_files = []

for reference in (references := scrnaseq_processed_data_dir.iterdir()):
    name = reference.stem
    counts = seurat_to_h5_dir / name / "RNA_counts.h5"

    for method in methods:
        for target_count in target_counts:

            out_file_coexpr = results_dir / f'coexpression_scrnaseq_references/{name}/coexpression_{method}_{target_count}.parquet' 
            out_file_pos_rate = results_dir / f'coexpression_scrnaseq_references/{name}/positivity_rate_{method}_{target_count}.parquet'

            out_files.extend([out_file_coexpr,out_file_pos_rate])

            rule:
                name: f'coexpression_scrnaseq_references/{name}/{method}_{target_count}'
                input:
                    counts=counts,
                output:
                    out_file_coexpr=out_file_coexpr,
                    out_file_pos_rate=out_file_pos_rate,
                params:
                    method=method,
                    target_count=target_count,
                    genes=genes,
                threads: 1
                resources:
                    mem='50GB',
                    runtime='30m',
                conda:
                    "general_cuda"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file_coexpr})"

                    python workflow/scripts/coexpression_sample.py \
                    -i {input.counts} \
                    -m {params.method} \
                    -c {params.target_count} \
                    --outcoexp {output.out_file_coexpr} \
                    --outposrate {output.out_file_pos_rate} \
                    -g {params.genes}

                    echo "DONE"
                    """


rule coexpression_scrnaseq_references_all:
    input:
        out_files
    output:
        touch(results_dir / "coexpression_scrnaseq_references_all.done")

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])
palette_dir = Path(config['xenium_metadata_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])
seurat_to_h5_dir = results_dir / 'seurat_to_h5'
scrnaseq_processed_data_dir = Path(config['scrnaseq_processed_data_dir'])

# Params
n_comps = config['umap_n_comps']
n_neighbors = config['umap_n_neighbors']
min_dist = config['umap_min_dist']
metric = config['umap_metric']

s=3
alpha=0.5
dpi = 300
points_only = True

cell_type_palette = palette_dir / 'col_palette_cell_types_combo.csv'
panel_palette = palette_dir / 'col_palette_panel.csv'
sample_palette = palette_dir / 'col_palette_sample.csv'

layer = 'RNA_counts'
colors = ['sample','Level2.1']#['Level1','Level2','Level3','Level4','panel','sample',] # condition and sample as color to plot added here in addition to levels
extension = 'png'

reference1 = scrnaseq_processed_data_dir / 'matched_combo_standard_lung_specific'
reference1_name = reference1.stem
reference1_dir = seurat_to_h5_dir / reference1_name

reference2 = scrnaseq_processed_data_dir / 'matched_combo_standard_breast_specific'
reference2_name = reference2.stem
reference2_dir = seurat_to_h5_dir / reference2_name

out_files_panel = []
for color in colors:

    # input embedding file (doesn't depend on ref,method or color loops but more readable to have here)
    embed_file = results_dir / f'embed_panel_scrnaseq_joint/{reference1_name}_{reference2_name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 

    # no need to plot panel for panel color UMAPs
    if color == 'panel':
        continue
    
    # no need to plot sample coloring for every param combination
    if color == 'sample':
        continue

    out_file = figures_dir / f"embed_panel_scrnaseq_joint/{reference1_name}_{reference2_name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}_{color}.{extension}"
    out_files_panel.append(out_file)

    rule:
        name: f'embed_panel_scrnaseq_joint_plot/{reference1_name}_{reference2_name}/umap_{layer}_{color}'
        input:
            embed_file=embed_file,
        output:
            out_file=out_file,
        params:
            reference1=reference1_dir,
            reference2=reference2_dir,
            color=color,
            cell_type_palette=cell_type_palette,
            panel_palette=panel_palette,
            sample_palette=sample_palette,
            s=s,
            alpha=alpha,
            dpi=dpi,
            points_only='--points_only' if points_only else '',
        threads: 1
        resources:
            mem='30GB',
            runtime='10m',
        conda:
            "general_cuda"
        shell:
            """
            mkdir -p "$(dirname {output.out_file})"

            python workflow/scripts/scRNAseq/embed_panel_scrnaseq_joint_plot.py \
            --embed_file {input.embed_file} \
            --reference1 {params.reference1} \
            --reference2 {params.reference2} \
            --color {params.color} \
            --out_file {output.out_file} \
            --cell_type_palette {params.cell_type_palette} \
            --panel_palette {params.panel_palette} \
            --sample_palette {params.sample_palette} \
            --s {params.s} \
            --alpha {params.alpha} \
            --dpi {params.dpi} \
            {params.points_only} \

            echo "DONE"
            """

rule embed_panel_scrnaseq_joint_plot_all:
    input:
        out_files_panel

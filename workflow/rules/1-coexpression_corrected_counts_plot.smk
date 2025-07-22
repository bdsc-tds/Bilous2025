
# take arbitrary segmentation just to loop over all conditions and panel combinations
segmentation = list(xenium_processed_data_dir.iterdir())[0] 

methods = ['conditional','jaccard']#,'pearson','spearman']
target_counts = [30,50,200]
min_positivity_rate = 0.01
min_cond_coex = 'auto'
min_cond_coex_mode = 'both'
cc_cutoff = 1.5
ref_segmentation = '10x_0um'
ref_oversegmentation = '10x_15um'
showfliers = False
log_scale = True
n_top_gene_pairs = 100_000_000 # big number to plot all gene pairs

out_files_panel = []

for correction_method in correction_methods:
    coexpression_corrected_counts_dir = results_dir / f'{correction_method}_coexpression'
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for method in methods:
                for target_count in target_counts:
                    if target_count > 50 and panel.stem != '5k':
                        continue

                    k = (condition.stem,panel.stem)
                    name = '/'.join(k)
                    min_cond_coex_str = f'{min_cond_coex=}' if isinstance(min_cond_coex,float) else min_cond_coex

                    panel_coexpression_corrected_counts = coexpression_corrected_counts_dir / name
                    out_file_plot_sample = figures_dir / f'{correction_method}_coexpression_panel/{name}/coexpression_{method}_{target_count=}_{min_cond_coex_str}_sample.{extension}'
                    out_file_plot_panel = figures_dir / f'{correction_method}_coexpression_panel/{name}/coexpression_{method}_{target_count=}_{min_cond_coex_str}_panel.{extension}'
                    out_file_gene_pairs = results_dir / f'{correction_method}_coexpression_gene_pairs/{name}/coexpression_gene_pairs_{method}_{target_count=}_{min_cond_coex_str}.parquet'
                    out_files_panel.extend([out_file_plot_sample,out_file_plot_panel,out_file_gene_pairs])

                    # adapt resources
                    if panel.stem == '5k':
                        if target_count > 50:
                            mem = '120GB'
                        else:
                            mem = '80GB'
                        runtime = '80m'
                    else:
                        mem = '20GB'
                        runtime = '10m'

                    rule:
                        name: f'{correction_method}_coexpression_plot_panel/{name}/coexpression_{method}_{target_count=}'
                        input:
                            coexpression_corrected_counts_is_done=results_dir / f"{correction_method}_coexpression.done",
                        output:
                            out_file_plot_sample=out_file_plot_sample,
                            out_file_plot_panel=out_file_plot_panel,
                            out_file_gene_pairs=out_file_gene_pairs,
                        params:
                            coexpression_corrected_counts_dir=coexpression_corrected_counts_dir,  
                            plot_condition=condition.stem,
                            plot_panel=panel.stem,
                            method=method,
                            target_count=target_count,
                            min_positivity_rate=min_positivity_rate,
                            min_cond_coex=min_cond_coex,
                            min_cond_coex_mode=min_cond_coex_mode,
                            cc_cutoff=cc_cutoff,
                            ref_segmentation=ref_segmentation,
                            ref_oversegmentation=ref_oversegmentation,
                            segmentation_palette=segmentation_palette,
                            dpi=dpi,
                            showfliers='--showfliers' if showfliers else '',
                            log_scale='--log_scale' if log_scale else '',
                        threads: 1
                        resources:
                            mem=mem,
                            runtime=runtime,
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file_plot_sample})"

                            python workflow/scripts/xenium/coexpression_panel_plot.py \
                            --coexpression_dir {params.coexpression_corrected_counts_dir} \
                            --plot_condition {params.plot_condition} \
                            --plot_panel {params.plot_panel} \
                            --out_file_plot_sample {output.out_file_plot_sample} \
                            --out_file_plot_panel {output.out_file_plot_panel} \
                            --out_file_gene_pairs {output.out_file_gene_pairs} \
                            --method {params.method} \
                            --target_count {params.target_count} \
                            --min_positivity_rate {params.min_positivity_rate} \
                            --cc_cutoff {params.cc_cutoff} \
                            --min_cond_coex {params.min_cond_coex} \
                            --min_cond_coex_mode {params.min_cond_coex_mode} \
                            --ref_segmentation {params.ref_segmentation} \
                            --ref_oversegmentation {params.ref_oversegmentation} \
                            --segmentation_palette {params.segmentation_palette} \
                            --dpi {params.dpi} \
                            {params.showfliers} \
                            {params.log_scale} \

                            echo "DONE"
                            """


out_files_conditions = []

for correction_method in correction_methods:
    for method in methods:
        for target_count in target_counts:
            if target_count > 50 and panel.stem != '5k':
                continue

            min_cond_coex_str = f'{min_cond_coex=}' if isinstance(min_cond_coex,float) else min_cond_coex

            out_file_plot = figures_dir / f'{correction_method}_coexpression_conditions/coexpression_{method}_{target_count=}_{min_cond_coex_str}_conditions.{extension}'
            out_file_gene_pairs = results_dir / f'{correction_method}_coexpression_conditions_gene_pairs/coexpression_gene_pairs_{method}_{target_count=}_{min_cond_coex_str}.parquet'
            out_files_conditions.extend([out_file_plot_sample,out_file_plot_panel,out_file_gene_pairs])

            # adapt resources
            if panel.stem == '5k':
                if target_count > 50:
                    mem = '120GB'
                else:
                    mem = '80GB'
                runtime = '80m'
            else:
                mem = '20GB'
                runtime = '10m'

            rule:
                name: f'{correction_method}_coexpression_conditions/coexpression_{method}_{target_count=}'
                input:
                    coexpression_corrected_counts_is_done=results_dir / f"{correction_method}_coexpression.done",
                output:
                    out_file_plot=out_file_plot,
                    out_file_gene_pairs=out_file_gene_pairs,
                params:
                    coexpression_corrected_counts_dir=coexpression_corrected_counts_dir,
                    method=method,
                    target_count=target_count,
                    min_positivity_rate=min_positivity_rate,
                    min_cond_coex=min_cond_coex,
                    min_cond_coex_mode=min_cond_coex_mode,
                    cc_cutoff=cc_cutoff,
                    ref_segmentation=ref_segmentation,
                    ref_oversegmentation=ref_oversegmentation,
                    segmentation_palette=segmentation_palette,
                    dpi=dpi,
                    showfliers='--showfliers' if showfliers else '',
                    log_scale='--log_scale' if log_scale else '',
                    n_top_gene_pairs=n_top_gene_pairs,
                threads: 1
                resources:
                    mem=mem,
                    runtime=runtime,
                conda:
                    "spatial"
                shell:
                    """
                    mkdir -p "$(dirname {output.out_file_plot})"

                    python workflow/scripts/xenium/coexpression_conditions_plot.py \
                    --coexpression_dir {params.coexpression_corrected_counts_dir} \
                    --out_file_plot {output.out_file_plot} \
                    --out_file_gene_pairs {output.out_file_gene_pairs} \
                    --method {params.method} \
                    --target_count {params.target_count} \
                    --min_positivity_rate {params.min_positivity_rate} \
                    --cc_cutoff {params.cc_cutoff} \
                    --min_cond_coex {params.min_cond_coex} \
                    --min_cond_coex_mode {params.min_cond_coex_mode} \
                    --ref_segmentation {params.ref_segmentation} \
                    --ref_oversegmentation {params.ref_oversegmentation} \
                    --segmentation_palette {params.segmentation_palette} \
                    --dpi {params.dpi} \
                    {params.showfliers} \
                    {params.log_scale} \
                    --n_top_gene_pairs {params.n_top_gene_pairs} \

                    echo "DONE"
                    """


rule coexpression_corrected_counts_plot_panel_all:
    input:
        out_files_panel

rule coexpression_corrected_counts_plot_conditions_all:
    input:
        out_files_conditions



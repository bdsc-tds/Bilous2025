import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.table import Table


def QC(sdata, donor_name="", save_path=None):
    adata = sdata["table"]
    metrics_summary = adata.uns["metrics_summary"].iloc[0]

    # precomputed QC metrics
    num_cells_detected = metrics_summary["num_cells_detected"]
    region_area = int(metrics_summary["region_area"])
    total_cell_area = int(metrics_summary["total_cell_area"])
    total_hq_transcripts = metrics_summary["total_high_quality_decoded_transcripts"]
    median_genes_per_cell = int(metrics_summary["median_genes_per_cell"])
    median_transcripts_per_cell = int(metrics_summary["median_transcripts_per_cell"])
    perc_transcripts_decoded_q20 = (
        metrics_summary["fraction_transcripts_decoded_q20"] * 100
    )
    estimated_number_of_false_positive_transcripts_per_cell = metrics_summary[
        "estimated_number_of_false_positive_transcripts_per_cell"
    ]

    # extra QC metrics on HQ genes
    adata.obs["n_genes_by_counts"] = (adata.X > 0).sum(1).A1
    hq_gene_transcripts = sdata["transcripts"][
        sdata["transcripts"]["feature_name"].isin(adata.var_names)
        & (sdata["transcripts"]["qv"] > 20)
    ]
    total_hq_gene_transcripts = len(hq_gene_transcripts.index)
    total_hq_transcripts_by_gene = (
        hq_gene_transcripts["feature_name"].value_counts().compute()
    )

    df_unassigned_hq_gene_transcripts = hq_gene_transcripts[
        hq_gene_transcripts["cell_id"] == "UNASSIGNED"
    ]

    unassigned_hq_gene_transcripts = len(df_unassigned_hq_gene_transcripts.index)
    unassigned_hq_transcripts_by_gene = (
        df_unassigned_hq_gene_transcripts["feature_name"].value_counts().compute()
    )
    perc_unassigned_hq_transcripts = (
        unassigned_hq_gene_transcripts / total_hq_gene_transcripts * 100
    )
    perc_unassigned_hq_transcripts_by_gene = (
        unassigned_hq_transcripts_by_gene / total_hq_transcripts_by_gene * 100
    )

    # probes QC
    perc_cprobes = (
        adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
    )
    perc_cwords = (
        adata.obs["control_codeword_counts"].sum()
        / adata.obs["total_counts"].sum()
        * 100
    )

    # Create a dictionary for title content to format into a table
    title_content = {
        "QC metrics": [
            "Number of cells",
            "Region area",
            "Total cell area",
            "Total high quality transcripts (qv>20)",
            "Median genes per cell",
            "Median transcripts per cell",
            "% high quality transcripts",
            "% unassigned high quality transcripts",
            "N° false positive transcripts per cell",
            "Negative DNA probe count %",
            "Negative decoding count %",
        ],
        "": [
            f"{num_cells_detected}",
            f"{region_area}",
            f"{total_cell_area}",
            f"{total_hq_transcripts}",
            f"{median_genes_per_cell}",
            f"{median_transcripts_per_cell:.2f}",
            f"{perc_transcripts_decoded_q20:.2f}",
            f"{perc_unassigned_hq_transcripts:.2f}",
            f"{estimated_number_of_false_positive_transcripts_per_cell:.3f}",
            f"{perc_cprobes:.3f}",
            f"{perc_cwords:.3f}",
        ],
    }

    # Create a figure and axis for the table
    fig, axs = plt.subplots(1, 5, figsize=(25, 6))
    axs = axs.flat

    # Add a table for the title
    fig.subplots_adjust(top=0.8)
    table_ax = fig.add_axes([0.32, 0.8, 0.4, 0.8])  # x, y, width, height
    table_ax.axis("off")

    # Create a Matplotlib table
    table = Table(table_ax, bbox=[0.2, 0.2, 0.55, 0.55])
    col_widths = [0.8, 0.2]  # Proportional widths for QC metric and Value columns

    # Add rows to the table
    for i, (metric, value) in enumerate(
        zip(title_content["QC metrics"], title_content[""])
    ):
        table.add_cell(
            i,
            0,
            col_widths[0],
            1 / len(title_content["QC metrics"]),
            text=metric,
            loc="left",
        )
        table.add_cell(
            i,
            1,
            col_widths[1],
            1 / len(title_content["QC metrics"]),
            text=value,
            loc="center",
        )

    # Add column headers
    table.add_cell(
        -1,
        0,
        col_widths[0],
        1 / len(title_content["QC metrics"]),
        text=donor_name,
        loc="center",
        facecolor="lightblue",
    )
    table.add_cell(
        -1,
        1,
        col_widths[1],
        1 / len(title_content["QC metrics"]),
        text="",
        loc="center",
        facecolor="lightblue",
    )
    # Adjust font size for table rows and headers
    for key, cell in table.get_celld().items():
        cell.set_fontsize(20)  # Adjust font size here

    # Attach table to axis
    table_ax.add_table(table)

    axs[0].set_title("Area of segmented cells")
    sns.histplot(
        adata.obs["cell_area"],
        kde=False,
        ax=axs[0],
    )

    axs[1].set_title("Nucleus ratio")
    sns.histplot(
        adata.obs["nucleus_area"] / adata.obs["cell_area"],
        kde=False,
        ax=axs[1],
    )
    axs[1].set_xlabel("nucleus ratio")

    axs[2].set_title("Total HQ transcripts per cell")
    sns.histplot(
        adata.obs["transcript_counts"],
        kde=False,
        ax=axs[2],
    )

    axs[3].set_title("Unique HQ transcripts per cell")
    sns.histplot(
        adata.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[3],
    )

    axs[4].set_title("% unassigned HQ transcripts by gene")
    sns.histplot(
        perc_unassigned_hq_transcripts_by_gene,
        kde=False,
        ax=axs[4],
    )
    axs[4].set_xlabel("percentage")

    for i, ax in enumerate(axs):
        if i == 0:
            ax.yaxis.get_label().set_fontsize(15)
        else:
            ax.set(ylabel="")
        ax.xaxis.get_label().set_fontsize(15)
        ax.title.set_fontsize(13)

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close()
    else:
        plt.show()

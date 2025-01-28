# convert_seurat_to_soma.R
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript convert_seurat_to_soma.R <seurat_rds> <soma_output>")
}

library(Seurat)
library(tiledbsoma)

seurat_rds <- args[1]
soma_output <- args[2]


seurat_rds <- '/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_0um/breast/breast/0OE1/0OE1/std_seurat_objects/preprocessed_seurat.rds'
soma_output <- 'workflow/scratch/test.soma'
seu <- readRDS(seurat_rds)
tiledbsoma::write_soma(seu, soma_output)

cat("Finished R tasks: Saved SOMA at", soma_output, "\n")

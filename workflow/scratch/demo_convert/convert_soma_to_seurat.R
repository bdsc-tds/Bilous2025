# convert_seurat_to_soma.R
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript convert_seurat_to_soma.R <soma_input> <seurat_output>")
}

library(Seurat)
library(tiledbsoma)

soma_input <- args[1]
seurat_output <- args[2]

experiment <- SOMAExperimentOpen(soma_input)
query <- SOMAExperimentAxisQuery$new(
  experiment = experiment,
  measurement_name = "RNA"
)

X_layers <- c()
for (name in query$ms$X$names()) {
    X_layers[name] <- name
}

seu <- query$to_seurat(
  X_layers=X_layers,
  obs_index = "obs_id",
  var_index = "var_id"
)

saveRDS(seu, seurat_output)

cat("Finished R tasks: Saved SOMA at", seurat_output, "\n")

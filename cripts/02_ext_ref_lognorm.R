# ============================================================
# Title: Log-normalization of Reference Objects
# Description: Applies Seurat LogNormalize to all references
# ============================================================

library(Seurat)

# ---------------------------
# Config
# ---------------------------
CONFIG <- list(
  ref_dir = "data/references",
  normalization_method = "LogNormalize"
)

options(future.globals.maxSize = 40000 * 1024^2)

# ---------------------------
# Files to process
# ---------------------------
files <- c(
  "external_lung.rds",
  "external_breast.rds"
)

# ---------------------------
# Helper: normalize object
# ---------------------------
normalize_seurat <- function(path) {

  obj <- readRDS(path)

  obj <- NormalizeData(
    obj,
    normalization.method = CONFIG$normalization_method,
    assay = "RNA"
  )

  return(obj)
}

# ---------------------------
# Run normalization
# ---------------------------
objects <- lapply(files, function(f) {
  in_path  <- file.path(CONFIG$ref_dir, f)
  out_path <- file.path(CONFIG$ref_dir, f)

  message("Processing: ", f)

  obj <- normalize_seurat(in_path)

  saveRDS(obj, out_path)

  return(NULL)
})

# ============================================================
# End
# ============================================================

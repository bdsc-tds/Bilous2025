library(Seurat)
library(dplyr)

# ── USER CONFIG ───────────────────────────────────────────────────────────────
LUNG_RDS    <- "/path/to/external_lung.rds"
BREAST_RDS  <- "/path/to/external_breast.rds"
LUNG_OUT    <- "/path/to/external_lung_sct.rds"
BREAST_OUT  <- "/path/to/external_breast_sct.rds"
LUNG_TARGET <- 500000   # cells to subsample for lung before SCTransform
SEED        <- 1
# ─────────────────────────────────────────────────────────────────────────────

options(future.globals.maxSize = 48000 * 1024^2)

subset_proportional <- function(seu, target_n, celltype_col = "Level3", seed = 1) {
  set.seed(seed)
  target_n  <- min(target_n, ncol(seu))
  meta      <- seu@meta.data
  cell_types <- meta[[celltype_col]]
  props     <- table(cell_types) / length(cell_types)
  sizes     <- round(props * target_n)
  cells     <- unlist(lapply(names(sizes), function(ct) {
    pool <- rownames(meta)[cell_types == ct]
    sample(pool, min(length(pool), sizes[[ct]]), replace = FALSE)
  }))
  subset(seu, cells = cells)
}

# ── Breast: SCTransform on full object ────────────────────────────────────────
breast <- readRDS(BREAST_RDS)
breast <- SCTransform(breast, verbose = FALSE)
saveRDS(breast, BREAST_OUT)
cat("Breast saved:", BREAST_OUT, "\n")
rm(breast); gc()

# ── Lung: proportional subsample then SCTransform ─────────────────────────────
lung <- readRDS(LUNG_RDS)
lung <- subset_proportional(lung, target_n = LUNG_TARGET, seed = SEED)
cat("Lung subsampled to", ncol(lung), "cells\n")
lung <- SCTransform(lung, verbose = FALSE)
saveRDS(lung, LUNG_OUT)
cat("Lung saved:", LUNG_OUT, "\n")

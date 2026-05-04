# ============================================================
# Title: Balanced Subsampling + SCTransform
# ============================================================

library(Seurat)
library(dplyr)

CONFIG <- list(
  ref_dir = "data/references"
)

# ---------------------------
# Function: proportional subsampling
# ---------------------------
subset_seurat_by_proportion <- function(seurat_obj,
                                       target_n,
                                       celltype_col = "cell_type",
                                       seed = 1) {

  set.seed(seed)

  total_cells <- ncol(seurat_obj)
  target_n <- min(target_n, total_cells)

  meta <- seurat_obj@meta.data
  cell_types <- meta[[celltype_col]]

  # proportions per cell type
  prop_table <- table(cell_types) / length(cell_types)

  # expected sample sizes
  sample_sizes <- round(prop_table * target_n)

  sampled_cells <- unlist(lapply(names(sample_sizes), function(ct) {

    cells_ct <- rownames(meta)[cell_types == ct]
    n_take <- min(length(cells_ct), sample_sizes[[ct]])

    sample(cells_ct, n_take, replace = FALSE)
  }))

  subset(seurat_obj, cells = sampled_cells)
}

# ---------------------------
# Load objects
# ---------------------------
breast <- readRDS(file.path(CONFIG$ref_dir,"external_breast.rds"))
seu    <- readRDS(file.path(CONFIG$ref_dir,"external_lung.rds"))

options(future.globals.maxSize = 48000 * 1024^2)

# ---------------------------
# SCTransform breast (full)
# ---------------------------
breast <- SCTransform(breast, verbose = FALSE)
saveRDS(breast, file.path(CONFIG$ref_dir,"external_breast.rds"))

# ---------------------------
# Subsample + SCTransform lung
# ---------------------------
seu_filt <- subset_seurat_by_proportion(
  seu,
  target_n = 500000,
  celltype_col = "Level3"
)

seu_filt <- SCTransform(seu_filt, verbose = FALSE)

# ---------------------------
# Checks
# ---------------------------
dim(seu_filt)
table(seu_filt$Level3)

# ---------------------------
# Save
# ---------------------------
saveRDS(seu_filt, file.path(CONFIG$ref_dir,"external_lung.rds")

# ============================================================
# End
# ============================================================

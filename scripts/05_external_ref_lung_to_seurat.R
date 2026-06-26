library(Seurat)
library(Matrix)
library(dplyr)

# ── USER CONFIG ───────────────────────────────────────────────────────────────
MTX_DIR <- "/path/to/external_lung_mtx"   # output directory from external_ref_lung_h5ad_to_mtx.py
OUT_RDS <- "/path/to/external_lung.rds"
# ─────────────────────────────────────────────────────────────────────────────

counts <- ReadMtx(
  mtx            = file.path(MTX_DIR, "matrix.mtx"),
  cells          = file.path(MTX_DIR, "barcodes.txt"),
  features       = file.path(MTX_DIR, "features.txt"),
  feature.column = 1
)

meta <- read.csv(file.path(MTX_DIR, "metadata.csv"), row.names = 1)

seu <- CreateSeuratObject(counts = counts, assay = "RNA")
new_cols <- setdiff(colnames(meta), colnames(seu@meta.data))
seu@meta.data <- cbind(seu@meta.data, meta[rownames(seu@meta.data), new_cols, drop = FALSE])

umap_path <- file.path(MTX_DIR, "umap.csv")
if (file.exists(umap_path)) {
  umap <- read.csv(umap_path, row.names = 1)
  seu[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(umap[colnames(seu), , drop = FALSE]),
    key        = "UMAP_",
    assay      = "RNA"
  )
}

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Level3: recode cell_type, consolidating endothelial and fibroblast subtypes
seu@meta.data <- seu@meta.data %>%
  mutate(Level3 = recode(cell_type,
    "conventional dendritic cell"     = "CD141-positive myeloid dendritic cell",
    "capillary endothelial cell"       = "endothelial cell",
    "vein endothelial cell"            = "endothelial cell",
    "pulmonary artery endothelial cell"= "endothelial cell",
    "bronchus fibroblast of lung"      = "fibroblast of lung",
    "malignant cell"                   = "malignant cell of lung"
  ))

# Level2: broad groupings from recoded Level3
seu@meta.data <- seu@meta.data %>%
  mutate(Level2 = recode(Level3,
    "alveolar macrophage"                      = "macrophage",
    "CD1c-positive myeloid dendritic cell"     = "dendritic cell",
    "CD141-positive myeloid dendritic cell"    = "dendritic cell",
    "CD4-positive, alpha-beta T cell"          = "T cell",
    "CD8-positive, alpha-beta T cell"          = "T cell",
    "classical monocyte"                       = "monocyte",
    "club cell"                                = "epithelial cell of lung",
    "dendritic cell"                           = "dendritic cell",
    "endothelial cell of lymphatic vessel"     = "endothelial cell",
    "mesothelial cell"                         = "mesothelial cell",
    "multiciliated epithelial cell"            = "epithelial cell of lung",
    "neutrophil"                               = "granulocyte",
    "non-classical monocyte"                   = "monocyte",
    "plasmacytoid dendritic cell"              = "dendritic cell",
    "pulmonary alveolar type 1 cell"           = "epithelial cell of lung",
    "pulmonary alveolar type 2 cell"           = "epithelial cell of lung",
    "regulatory T cell"                        = "T cell"
  ))

# Level1: coarsest groupings from Level2
seu@meta.data <- seu@meta.data %>%
  mutate(Level1 = recode(Level2,
    "dendritic cell"        = "myeloid cell",
    "endothelial cell"      = "stromal cell",
    "epithelial cell of lung" = "epithelial cell of lung",
    "fibroblast of lung"    = "stromal cell",
    "granulocyte"           = "myeloid cell",
    "macrophage"            = "myeloid cell",
    "mesothelial cell"      = "stromal cell",
    "monocyte"              = "myeloid cell",
    "myeloid cell"          = "myeloid cell",
    "pericyte"              = "stromal cell",
    "plasma cell"           = "B cell",
    "smooth muscle cell"    = "stromal cell",
    "stromal cell"          = "stromal cell"
  ))

cat("Level1:\n"); print(table(seu$Level1))
cat("Level2:\n"); print(table(seu$Level2))
cat("Level3:\n"); print(table(seu$Level3))

saveRDS(seu, OUT_RDS)
cat("Saved:", OUT_RDS, "\n")
cat("Cells:", ncol(seu), "| Genes:", nrow(seu), "\n")

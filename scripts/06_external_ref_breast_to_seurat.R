library(Seurat)
library(Matrix)
library(dplyr)
library(readr)

# ── USER CONFIG ───────────────────────────────────────────────────────────────
# Source: Wu et al. 2021 BRCA scRNA-seq (barcodes.tsv / genes.tsv / matrix.mtx + metadata.csv)
COUNTS_DIR   <- "/path/to/Wu_etal_2021_BRCA_scRNASeq"   # 10X directory
METADATA_CSV <- "/path/to/Wu_etal_2021_BRCA_scRNASeq/metadata.csv"
OUT_RDS      <- "/path/to/external_breast.rds"
# ─────────────────────────────────────────────────────────────────────────────

counts <- Read10X(COUNTS_DIR, gene.column = 1)
meta   <- read.csv(METADATA_CSV, row.names = 1)

seu <- CreateSeuratObject(counts = counts, assay = "RNA")
new_cols <- setdiff(colnames(meta), colnames(seu@meta.data))
seu@meta.data <- cbind(seu@meta.data, meta[rownames(seu@meta.data), new_cols, drop = FALSE])

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Level3: recode celltype_minor
seu@meta.data <- seu@meta.data %>%
  mutate(Level3 = recode(celltype_minor,
    "Endothelial ACKR1"                              = "endothelial cell",
    "Endothelial RGS5"                               = "endothelial cell",
    "Endothelial CXCL12"                             = "endothelial cell",
    "Endothelial Lymphatic LYVE1"                    = "endothelial cell of lymphatic vessel",
    "CAFs MSC iCAF-like"                             = "fibroblast of breast",
    "CAFs myCAF-like"                                = "fibroblast of breast",
    "PVL Differentiated"                             = "pericyte",
    "PVL Immature"                                   = "pericyte",
    "Cycling PVL"                                    = "pericyte",
    "B cells Memory"                                 = "B cell",
    "B cells Naive"                                  = "B cell",
    "T cells CD8+"                                   = "CD8-positive, alpha-beta T cell",
    "T cells CD4+"                                   = "CD4-positive, alpha-beta T cell",
    "NK cells"                                       = "natural killer cell",
    "Cycling T-cells"                                = "cycling T cell",
    "NKT cells"                                      = "natural killer T cell",
    "Macrophage"                                     = "macrophage",
    "Monocyte"                                       = "monocyte",
    "Cycling_Myeloid"                                = "cycling myeloid cell",
    "DCs"                                            = "dendritic cell",
    "Myoepithelial"                                  = "myoepithelial cell",
    "Luminal Progenitors"                            = "progenitor cell of mammary luminal epithelium",
    "Mature Luminal"                                 = "luminal epithelial cell of mammary gland",
    "Plasmablasts"                                   = "plasma cell",
    "Cancer Cycling"                                 = "cycling malignant cell",
    "Cancer Her2 SC"                                 = "malignant cell HER2",
    "Cancer LumB SC"                                 = "malignant cell LUMB",
    "Cancer Basal SC"                                = "malignant cell basal",
    "Cancer LumA SC"                                 = "malignant cell LUMA"
  ))

# Level2: broad groupings from Level3
seu@meta.data <- seu@meta.data %>%
  mutate(Level2 = recode(Level3,
    "endothelial cell of lymphatic vessel"                = "endothelial cell",
    "CD8-positive, alpha-beta T cell"                    = "T cell",
    "CD4-positive, alpha-beta T cell"                    = "T cell",
    "natural killer T cell"                              = "T cell",
    "cycling T cell"                                     = "T cell",
    "cycling myeloid cell"                               = "myeloid cell",
    "myoepithelial cell"                                 = "epithelial cell of breast",
    "progenitor cell of mammary luminal epithelium"      = "epithelial cell of breast",
    "luminal epithelial cell of mammary gland"           = "epithelial cell of breast",
    "cycling malignant cell"                             = "malignant cell of breast",
    "malignant cell HER2"                                = "malignant cell of breast",
    "malignant cell LUMB"                                = "malignant cell of breast",
    "malignant cell basal"                               = "malignant cell of breast",
    "malignant cell LUMA"                                = "malignant cell of breast"
  ))

# Level1: coarsest groupings from Level2
seu@meta.data <- seu@meta.data %>%
  mutate(Level1 = recode(Level2,
    "epithelial cell of breast" = "epithelial cell of breast",
    "endothelial cell"          = "stromal cell",
    "fibroblast of breast"      = "stromal cell",
    "pericyte"                  = "stromal cell",
    "macrophage"                = "myeloid cell",
    "monocyte"                  = "myeloid cell",
    "myeloid cell"              = "myeloid cell",
    "dendritic cell"            = "myeloid cell",
    "plasma cell"               = "B cell"
  ))

cat("Level1:\n"); print(table(seu$Level1))
cat("Level2:\n"); print(table(seu$Level2))
cat("Level3:\n"); print(table(seu$Level3))

saveRDS(seu, OUT_RDS)
cat("Saved:", OUT_RDS, "\n")
cat("Cells:", ncol(seu), "| Genes:", nrow(seu), "\n")

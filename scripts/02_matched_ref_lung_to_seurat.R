library(Seurat)
library(Matrix)
library(dplyr)

# ── USER CONFIG ───────────────────────────────────────────────────────────────
MTX_DIR <- "/path/to/lung_mtx"       # output directory from 01_h5ad_to_mtx.py
OUT_RDS <- "/path/to/lung.rds"
# ─────────────────────────────────────────────────────────────────────────────

mtx_dir <- MTX_DIR

counts <- ReadMtx(
  mtx            = file.path(mtx_dir, "matrix.mtx"),
  cells          = file.path(mtx_dir, "barcodes.txt"),
  features       = file.path(mtx_dir, "features.txt"),
  feature.column = 1
)

seu <- CreateSeuratObject(counts = counts, assay = "RNA")

meta <- read.csv(file.path(mtx_dir, "metadata.csv"), row.names = 1)
shared <- intersect(rownames(seu@meta.data), rownames(meta))
new_cols <- setdiff(colnames(meta), colnames(seu@meta.data))
seu@meta.data <- cbind(seu@meta.data, meta[shared, new_cols, drop = FALSE])

umap_path <- file.path(mtx_dir, "umap.csv")
if (file.exists(umap_path)) {
  umap <- read.csv(umap_path, row.names = 1)
  seu[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(umap[colnames(seu), , drop = FALSE]),
    key        = "UMAP_",
    assay      = "RNA"
  )
}

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Level1 recoded from original Level3
seu@meta.data <- seu@meta.data %>%
  mutate(Level1 = recode(Level3,
    Epi_lung             = "epithelial cell of lung",
    Mast_cell            = "myeloid cell",
    Muscle_smooth        = "stromal cell",
    Macrophage           = "myeloid cell",
    B_cell               = "B cell",
    NK                   = "natural killer cell",
    T_CD4                = "T cell",
    T_CTL                = "T cell",
    Pericyte             = "stromal cell",
    Neutrophil           = "myeloid cell",
    DC_pc                = "myeloid cell",
    B_plasma             = "B cell",
    T_reg                = "T cell",
    Monocyte             = "myeloid cell",
    DC_1                 = "myeloid cell",
    Tu_B4                = "malignant cell of breast",
    Tu_B3                = "malignant cell of breast",
    Tu_B1                = "malignant cell of breast",
    Tu_L4                = "malignant cell of lung",
    Tu_L2                = "malignant cell of lung",
    Tu_L3                = "malignant cell of lung",
    Tu_L1                = "malignant cell of lung",
    Endothelia_vascular  = "stromal cell",
    DC_2                 = "myeloid cell",
    Fibroblast           = "stromal cell",
    Fibroblast_1GVR      = "stromal cell",
    Granulocyte          = "myeloid cell",
    T_CD8_exhausted      = "T cell",
    T_CXCL13            = "T cell",
    TNK_dividing         = "cycling lymphocyte",
    DC_activated         = "myeloid cell"
  ))

# Level2 recoded from original Level3
seu@meta.data <- seu@meta.data %>%
  mutate(Level2 = recode(Level3,
    Epi_lung             = "epithelial cell of lung",
    Mast_cell            = "mast cell",
    Muscle_smooth        = "smooth muscle cell",
    Macrophage           = "macrophage",
    B_cell               = "B cell",
    NK                   = "natural killer cell",
    T_CD4                = "T cell",
    T_CTL                = "T cell",
    Pericyte             = "pericyte",
    Neutrophil           = "granulocyte",
    DC_pc                = "dendritic cell",
    B_plasma             = "plasma cell",
    T_reg                = "T cell",
    Monocyte             = "monocyte",
    DC_1                 = "dendritic cell",
    Tu_B4                = "malignant cell of breast",
    Tu_B3                = "malignant cell of breast",
    Tu_B1                = "malignant cell of breast",
    Tu_L4                = "malignant cell of lung",
    Tu_L2                = "malignant cell of lung",
    Tu_L3                = "malignant cell of lung",
    Tu_L1                = "malignant cell of lung",
    Endothelia_vascular  = "endothelial cell",
    DC_2                 = "dendritic cell",
    Fibroblast           = "fibroblast of lung",
    Fibroblast_1GVR      = "fibroblast of lung",
    Granulocyte          = "granulocyte",
    T_CD8_exhausted      = "T cell",
    T_CXCL13            = "T cell",
    TNK_dividing         = "cycling lymphocyte",
    DC_activated         = "dendritic cell"
  ))

# Level3 recoded in-place
seu@meta.data <- seu@meta.data %>%
  mutate(Level3 = recode(Level3,
    Epi_lung             = "epithelial cell of lung",
    Mast_cell            = "mast cell",
    Muscle_smooth        = "smooth muscle cell",
    Macrophage           = "macrophage",
    B_cell               = "B cell",
    NK                   = "natural killer cell",
    T_CD4                = "CD4-positive, alpha-beta T cell",
    T_CTL                = "CD8-positive, alpha-beta T cell",
    Pericyte             = "pericyte",
    Neutrophil           = "neutrophil",
    DC_pc                = "plasmacytoid dendritic cell",
    B_plasma             = "plasma cell",
    T_reg                = "regulatory T cell",
    Monocyte             = "monocyte",
    DC_1                 = "CD141-positive myeloid dendritic cell",
    Tu_B4                = "malignant cell B4",
    Tu_B3                = "malignant cell B3",
    Tu_B1                = "malignant cell B1",
    Tu_L4                = "malignant cell L4",
    Tu_L2                = "malignant cell L2",
    Tu_L3                = "malignant cell L3",
    Tu_L1                = "malignant cell L1",
    Endothelia_vascular  = "endothelial cell",
    DC_2                 = "CD1c-positive myeloid dendritic cell",
    Fibroblast           = "fibroblast of lung",
    Fibroblast_B3        = "fibroblast of breast",
    Granulocyte          = "granulocyte",
    T_CD8_exhausted      = "exhausted T cell",
    T_CXCL13            = "T cell CXCL13",
    TNK_dividing         = "cycling lymphocyte",
    DC_activated         = "activated dendritic cell"
  ))

# Level4: rename from Harmonised_Level4 then recode
seu@meta.data <- seu@meta.data %>%
  dplyr::rename(Level4 = "Harmonised_Level4") %>%
  mutate(Level4 = recode(Level4,
    Epi_lung                           = "epithelial cell of lung",
    Mast_cell                          = "mast cell",
    Muscle_smooth                      = "smooth muscle cell",
    Macrophage                         = "macrophage",
    B_cell                             = "B cell",
    NK                                 = "natural killer cell",
    T_CD4                              = "CD4-positive, alpha-beta T cell",
    T_CTL                              = "CD8-positive, alpha-beta T cell",
    Pericyte                           = "pericyte",
    Neutrophil                         = "neutrophil",
    DC_pc                              = "plasmacytoid dendritic cell",
    B_plasma                           = "plasma cell",
    T_reg                              = "regulatory T cell",
    Monocyte                           = "monocyte",
    DC_1                               = "CD141-positive myeloid dendritic cell",
    Endothelia_vascular                = "endothelial cell",
    DC_2                               = "CD1c-positive myeloid dendritic cell",
    Fibroblast                         = "fibroblast of lung",
    Fibroblast_1GVR                    = "fibroblast of lung",
    Granulocyte                        = "granulocyte",
    T_CD8_exhausted                    = "exhausted CD8-positive, alpha-beta T cell",
    T_CXCL13                          = "T cell CXCL13",
    TNK_dividing                       = "cycling lymphocyte",
    DC_activated                       = "activated dendritic cell",
    B_plasma_IGHA1                     = "IgA1 plasma cell",
    B_plasma_IGHG1                     = "IgG1 plasma cell",
    B_plasma_IGHG3                     = "IgG3 plasma cell",
    B_plasma_IGHM                      = "IgM plasma cell",
    B_plasma_IGKC                      = "IgKC plasma cell",
    B_plasma_IGLC1                     = "IgLC1 plasma cell",
    Tu_L1_SFTPB                        = "malignant cell L1 SFTPB",
    Tu_L2_FXYD2                        = "malignant cell L2 FXYD2",
    Tu_B4_RHOB                         = "malignant cell B4 RHOB",
    Tu_L3_G0S2                         = "malignant cell L3 G0S2",
    Tu_L4_KRT17_immune_signature       = "malignant cell L4 KRT17 immune signature",
    Tu_L4_KRT17_mucous                 = "malignant cell L4 KRT17 mucous",
    Tu_L4_KRT17_necrosis               = "malignant cell L4 KRT17 MT",
    Tu_L4_KRT17_neutrophil_signature   = "malignant cell L4 KRT17 neutrophil signature",
    Tu_B3_CYP4F8                       = "malignant cell B3 CYP4F8",
    Tu_L3_G0S2_immune_signature        = "malignant cell L3 G0S2 immune signatue",
    Tu_B1_MUCL1                        = "malignant cell B1 MUCL1",
    Tu_B1_MUCL1_necrosis               = "malignant cell B1 MUCL1 MT",
    Tu_B1_MUCL1_transcription          = "malignant cell B1 MUCL1 transcription"
  ))

cat("Level1:\n"); print(table(seu$Level1))
cat("Level2:\n"); print(table(seu$Level2))
cat("Level3:\n"); print(table(seu$Level3))
cat("Level4:\n"); print(table(seu$Level4))

out_rds <- OUT_RDS
saveRDS(seu, out_rds)
cat("Saved:", out_rds, "\n")
cat("Cells:", ncol(seu), "| Genes:", nrow(seu), "\n")

# ============================================================
# Title: Reference Preparation (Lung + Breast scRNA-seq)
# ============================================================

library(Seurat)
library(dplyr)
library(biomaRt)
library(readr)

# ---------------------------
# Config
# ---------------------------
CONFIG <- list(
  ref_dir = "data/references",
  raw_dir = "data/raw",
  mapping_file = "simplified_cell_type_mapping.csv"
)

# ---------------------------
# Helper: gene symbol mapping
# ---------------------------
map_ensembl_to_symbol <- function(counts_matrix) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = rownames(counts_matrix),
    mart = ensembl
  ) %>%
    filter(hgnc_symbol != "") %>%
    distinct(hgnc_symbol, .keep_all = TRUE)

  counts_filt <- counts_matrix[gene_info$ensembl_gene_id, ]
  rownames(counts_filt) <- gene_info$hgnc_symbol

  counts_filt
}

# ============================================================
# LUNG
# ============================================================

lung_ref <- readRDS(file.path(CONFIG$ref_dir, "lung_reference.rds"))
Idents(lung_ref) <- "cell_type"

counts <- lung_ref@assays$RNA@counts
counts_filt <- map_ensembl_to_symbol(counts)

# Load mapping table
level_mapping <- read_csv(CONFIG$mapping_file)

meta <- lung_ref@meta.data
colnames(meta)[colnames(meta) == "cell_type"] <- "Level3"

meta <- meta %>%
  left_join(level_mapping, by = c("Level3" = "original")) %>%
  rename(Level1 = level1, Level2 = level2)

# ---- KEEP ALL YOUR CURATION ----

meta <- meta %>%
  mutate(
    Level2 = recode(Level2,
      "epithelial cell" = "epithelial cell of lung",
      "neutrophil" = "granulocyte",
      "fibroblast" = "fibroblast of lung"
    ),
    Level3 = recode(Level3,
      "conventional dendritic cell" = "CD141-positive myeloid dendritic cell",
      "capillary endothelial cell" = "endothelial cell",
      "vein endothelial cell" = "endothelial cell",
      "pulmonary artery endothelial cell" = "endothelial cell",
      "bronchus fibroblast of lung" = "fibroblast of lung"
    )
  )

meta <- meta %>%
  mutate(across(c(Level1, Level2, Level3),
                ~ recode(., "malignant cell" = "malignant cell of lung")))

meta <- meta %>%
  mutate(
    Level1 = recode(Level2,
      "monocyte" = "myeloid cell",
      "macrophage" = "myeloid cell",
      "endothelial cell" = "stromal cell",
      "dendritic cell" = "myeloid cell",
      "epithelial cell of lung" = "epithelial cell",
      "fibroblast of lung" = "stromal cell",
      "plasma cell" = "B cell",
      "smooth muscle cell" = "stromal cell",
      "pericyte" = "stromal cell",
      "granulocyte" = "myeloid cell"
    )
  )

lung_seu <- CreateSeuratObject(counts_filt, meta.data = meta)

# ============================================================
# BREAST
# ============================================================

counts <- Read10X(file.path(CONFIG$raw_dir, "breast_counts"))
meta <- read_csv(file.path(CONFIG$raw_dir, "breast_metadata.csv"))

breast <- CreateSeuratObject(counts = counts, meta.data = meta)

# ---- FULL ORIGINAL MAPPINGS PRESERVED ----

breast@meta.data <- breast@meta.data %>%
  mutate(
    Level3 = recode(celltype_minor,
      "Endothelial ACKR1"= "endothelial cell",
      "Endothelial RGS5"= "endothelial cell",
      "Endothelial CXCL12"= "endothelial cell",
      "CAFs MSC iCAF-like"="fibroblast of breast",
      "CAFs myCAF-like"="fibroblast of breast",
      "PVL Differentiated"="pericyte",
      "PVL Immature"="pericyte",
      "Endothelial Lymphatic LYVE1"="endothelial cell of lymphatic vessel",
      "B cells Memory"="B cell",
      "B cells Naive"="B cell",
      "T cells CD8+"= "CD8-positive, alpha-beta T cell",
      "T cells CD4+"="CD4-positive, alpha-beta T cell",
      "NK cells"="natural killer cell",
      "Cycling T-cells"="cycling T cell",
      "NKT cells"= "natural killer T cell",
      "Macrophage"="macrophage",
      "Monocyte"="monocyte",
      "Cycling_Myeloid"="cycling myeloid cell",
      "DCs"="dendritic cell",
      "Myoepithelial"="myoepithelial cell",
      "Luminal Progenitors"="progenitor cell of mammary luminal epithelium",
      "Mature Luminal"="luminal epithelial cell of mammary gland",
      "Plasmablasts"="plasma cell",
      "Cancer Cycling"="cycling malignant cell",
      "Cancer Her2 SC"="malignant cell HER2",
      "Cancer LumB SC"="malignant cell LUMB",
      "Cancer Basal SC"="malignant cell basal",
      "Cycling PVL"="pericyte",
      "Cancer LumA SC"="malignant cell LUMA"
    )
  )

breast@meta.data <- breast@meta.data %>%
  mutate(
    Level2 = recode(Level3,
      "endothelial cell of lymphatic vessel"="endothelial cell",
      "CD8-positive, alpha-beta T cell"="T cell",
      "CD4-positive, alpha-beta T cell"="T cell",
      "natural killer T cell"="T cell",
      "cycling T cell"="T cell",
      "cycling myeloid cell"="myeloid cell",
      "myoepithelial cell"="epithelial cell of breast",
      "progenitor cell of mammary luminal epithelium"="epithelial cell of breast",
      "luminal epithelial cell of mammary gland"="epithelial cell of breast",
      "cycling malignant cell"="malignant cell of breast",
      "malignant cell HER2"="malignant cell of breast",
      "malignant cell LUMB"="malignant cell of breast",
      "malignant cell basal"="malignant cell of breast",
      "malignant cell LUMA"="malignant cell of breast"
    )
  )

breast@meta.data <- breast@meta.data %>%
  mutate(
    Level1 = recode(Level2,
      "epithelial cell of breast"="epithelial cell",
      "endothelial cell"="stromal cell",
      "fibroblast of breast"="stromal cell",
      "pericyte"="stromal cell",
      "macrophage"="myeloid cell",
      "monocyte"="myeloid cell",
      "dendritic cell"="myeloid cell",
      "plasma cell"="B cell"
    )
  )

# ============================================================
# SAVE
# ============================================================

saveRDS(breast, file.path(CONFIG$ref_dir, "external_breast.rds"))
saveRDS(lung_seu, file.path(CONFIG$ref_dir, "external_lung.rds"))

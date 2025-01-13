# Load required libraries
library(Seurat)
library(tidyverse)
library(scuttle)
library(biomaRt)

### Chrom for checks ###
chrom <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/xenium_paper/scRNAseq/processed/matched_ref_breast_lung.rds")

### LUNG ###
# Load and filter NSCLC reference data
ref <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/mesothelioma/scRNAseq/processed/references/cellxgene/NSCLC_atlas_extended.rds")
Idents(ref) <- "cell_type"
counts <- ref@assays$RNA@counts

# Fetch gene symbols using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(counts),
  mart = ensembl
)

gene_info_filt <- gene_info %>%
  filter(if_all(everything(), ~ . != ""))


# Filter the data frame to retain only rows with valid genes
valid_genes_df <- gene_info_filt[!is.na(gene_info_filt$hgnc_symbol),]
unique_genes <- valid_genes_df[!duplicated(valid_genes_df$hgnc_symbol),]

counts_filt <- counts[unique_genes[["ensembl_gene_id"]],]

rownames(counts_filt) <- unique_genes[["hgnc_symbol"]]

cat("Old reference dimensions:", dim(counts), "\n")
cat("New reference dimensions:", dim(counts_filt), "\n")

# Read the second CSV file for additional mappings
level_mapping <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/mesothelioma/results/annotation/simplified_cell_type_mapping.csv",
                          stringsAsFactors = FALSE)

# change cell_type to level3
colnames(ref@meta.data)[colnames(ref@meta.data) == "cell_type"] <- "Level3"

# Merge the main data frame with the level mapping table
ref@meta.data <- ref@meta.data %>%
  left_join(level_mapping, by = c("Level3" = "original"))
table(ref@meta.data$level1, useNA = "always")

colnames(ref@meta.data)[colnames(ref@meta.data) == "level1"] <- "Level1"
colnames(ref@meta.data)[colnames(ref@meta.data) == "level2"] <- "Level2"

ref@meta.data <- ref@meta.data %>%
  mutate(Level2 = recode(Level2,
                         "epithelial cell" = "epithelial cell of lung",
                         "mast cell" = "granulocyte",
                         "neutrophil" = "granulocyte",
                         "fibroblast" = "fibroblast of lung"
                         ))

ref@meta.data <- ref@meta.data %>%
  mutate(Level3 = recode(Level3,
                         "conventional dendritic cell" = "CD141-positive myeloid dendritic cell",
                         "capillary endothelial cell" = "endothelial cell",
                         "vein endothelial cell"= "endothelial cell",
                         "pulmonary artery endothelial cell"="endothelial cell"
                         ))

ref@meta.data <- ref@meta.data %>%
  mutate(across(c(Level1, Level2, Level3),
                ~ recode(.,
                         "malignant cell" = "malignant cell of lung")))

ref@meta.data <- ref@meta.data %>%
  mutate(Level1 = recode(Level2,
                         "monocyte" = "myeloid cell",
                         "macrophage" = "myeloid cell",             
                         "endothelial cell"  = "stromal cell",   
                         "natural killer cell"  = "T cell",   
                         "dendritic cell" = "myeloid cell",            
                         "epithelial cell of lung" = "epithelial cell",
                         "fibroblast of lung"= "stromal cell",      
                         "plasma cell" = "B cell",
                         "smooth muscle cell"= "stromal cell",               
                         "pericyte" = "stromal cell",              
                         "granulocyte"= "myeloid cell"
                         
  ))

#re-make seurat obj
seu_filt <- CreateSeuratObject(counts_filt, meta.data = ref@meta.data)
rownames(seu_filt@meta.data) 


### BREAST ### 
# make seurat object
path <- Read10X("/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/xenium_paper/scRNAseq/processed/Wu_etal_2021_BRCA_scRNASeq",
                gene.column = 1)
metadata <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/xenium_paper/scRNAseq/processed/Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
breast <- CreateSeuratObject(counts = path,
                             meta.data = metadata)

# standardise annotations
unique(breast$celltype_minor)

breast@meta.data <- breast@meta.data %>%
  mutate(Level3 = recode(celltype_minor,
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
                         "Cycling T-cells"="dividing TNK cell",
                         "NKT cells"= "natural killer T cell",
                         "Macrophage"="macrophage",
                         "Monocyte"="monocyte",
                         "Cycling_Myeloid"="dividing myeloid cell",
                         "DCs"="dendritic cell",
                         "Myoepithelial"="myoepithelial cell",
                         "Luminal Progenitors"="progenitor cell of mammary luminal epithelium",
                         "Mature Luminal"="luminal epithelial cell of mammary gland",
                         "Plasmablasts"="plasma cell",
                         "Cancer Cycling"="dividing malignant cell",
                         "Cancer Her2 SC"="malignant cell HER2",
                         "Cancer LumB SC"="malignant cell LUMB",
                         "Cancer Basal SC"="malignant cell basal",
                         "Cycling PVL"="pericyte",
                         "Cancer LumA SC"="malignant cell LUMA"
  ))
table(breast$Level3)


breast@meta.data <- breast@meta.data %>%
  mutate(Level2 = recode(Level3,
                         "endothelial cell of lymphatic vessel"="endothelial cell",
                         "B cell"="B cell" ,
                         "CD8-positive, alpha-beta T cell"="T cell",
                         "CD4-positive, alpha-beta T cell"="T cell",
                         "natural killer cell"="natural killer cell",
                         "dividing TNK cell"="natural killer cell",
                         "natural killer T cell"="T cell",
                         "dividing myeloid cell"="myeloid cell",
                         "dendritic cell"="dendritic cell",
                         "myoepithelial cell"="epithelial cell of breast",
                         "progenitor cell of mammary luminal epithelium"="epithelial cell of breast",
                         "luminal epithelial cell of mammary gland"="epithelial cell of breast",
                         "dividing malignant cell"="malignant cell of breast",
                         "malignant cell HER2"="malignant cell of breast",
                         "malignant cell LUMB"="malignant cell of breast",
                         "malignant cell basal"="malignant cell of breast",
                         "dividing pericyte"="pericyte",
                         "malignant cell LUMA"="malignant cell of breast"
  ))


breast@meta.data <- breast@meta.data %>%
  mutate(Level1 = recode(Level2,
                         "epithelial cell of breast"="epithelial cell",
                         "endothelial cell"="stromal cell",
                         "fibroblast of breast"="stromal cell",
                         "pericyte"="stromal cell",
                         "natural killer cell"="T cell",
                         "macrophage"="myeloid cell",
                         "monocyte"="myeloid cell",
                         "dendritic cell"="myeloid cell",
                         "plasma cell"="B cell"
  ))

### CHECKS ###
setdiff(chrom$Level1, breast$Level1)
setdiff(breast$Level1, chrom$Level1)

setdiff(chrom$Level1, seu_filt$Level1)
setdiff(seu_filt$Level1, chrom$Level1)


setdiff(chrom$Level2, breast$Level2)
setdiff(breast$Level2, chrom$Level2)

setdiff(chrom$Level2, seu_filt$Level2)
setdiff(seu_filt$Level2, chrom$Level2)


setdiff(chrom$Level3, breast$Level3)
setdiff(breast$Level3, chrom$Level3)

setdiff(chrom$Level3, seu_filt$Level3)
setdiff(seu_filt$Level3, chrom$Level3)

# Save the updated reference objects
saveRDS(breast, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/xenium_paper/scRNAseq/processed/breast_external.rds")
saveRDS(seu_filt, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/xenium_paper/scRNAseq/processed/lung_external.rds")


library(dplyr)
library(Seurat)
library(arrow)

script_dir <- this.path::here()
source(paste0(script_dir,"/../scripts/xenium/celladmix_decontam.R"))
source(paste0(script_dir,"/../scripts/xenium/celladmix_utils.R"))
source(paste0(script_dir,"/../scripts/xenium/celladmix_plots.R"))
source(paste0(script_dir,"/../scripts/xenium/celladmix_data_loading.R"))

TRANSCRIPT_PATH <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_0um/breast/breast/0OE1/0OE1/normalised_results/outs/transcripts.parquet"
RCTD_PATH <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_0um/breast/breast/0OE1/0OE1/cell_type_annotation/reference_based/matched_reference/rctd/Level1/single_cell/labels.csv"

# load the molecule-level spatial data
mols.qv.threshold <- 20
transcripts <- read_parquet(TRANSCRIPT_PATH)
transcripts <- filter(transcripts, qv >= mols.qv.threshold)
transcripts <- filter(transcripts, is_gene) # remove dummy molecules

# load the RCTD object to get cell type annotations
cell_meta <- read.csv(RCTD_PATH,row.names = 1)
cell_meta$cell <- rownames(cell_meta)

transcripts <- prepare_transcript_data_xenium(transcripts,cell_meta)
counts <- transcripts_to_count_matrix(transcripts)

# reduce data to just the fibroblast cells in the tested regions
ct <- 'B cell'
cells_keep <- cell_meta %>% filter(first_type==ct) %>% rownames()
cell_meta_keep <- cell_meta[cells_keep,]
counts_keep <- counts[,cells_keep]
transcripts_keep <- transcripts %>% filter(celltype==ct)

# subset to only fibroblasts inside the tumor interior niche to help pull out the contamination signal

match_idx <- match(transcripts_keep$cell,rownames(cell_meta))
# df_subs$niche <- cell_meta$niche[match_idx]
# df_tumor <- df_subs[df_subs$niche=='tumor interior',]

# removing unneeded columns from the df
# df_tumor <- df_tumor[,c('cell','gene','niche','mol_id','x','y','z')]

# compute the molecule NCV matrix
X <- get_knn_counts(transcripts_keep, k=20, ncores = 20) %>% as.matrix()
res <- run_weighted_nmf(X, k=5, n.downdonor=100)
expr <- res@fit@H %>% {. / rowSums(.)} %>% t() %>% magrittr::set_colnames(1:ncol(.))

malignant_admix_factor <- 4
# Next, we will run CRF to assign each molecule to one of these factors. Then, we can remove the molecules coming from contaminating cell types.

# run crf
df_fib_crf <- run_crf_per_cell(transcripts_keep, expr, ncores=30)

# this appended a column 'factor' indicating which factor each molecule has been assigned to
head(df_fib_crf, 3)

# now, remove the molecules assigned to the admixture factor
df_fib_cln <- filter(df_fib_crf, !(factor %in% malignant_admix_factor))

print(paste0('Number of molecules originally: ', nrow(df_fib_crf)))
print(paste0('Number of molecules after cleaning: ', nrow(df_fib_cln)))

# compute a clean counts matrix from molecule level data
dat_cln <- df_fib_cln %>% transcripts_to_count_matrix() %>% .[,rownames(cell_meta_keep)]

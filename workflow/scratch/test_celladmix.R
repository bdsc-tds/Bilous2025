library(dplyr)

script_dir <- this.path::here()
source(paste0(script_dir,"/../scripts/xenium/celladmix_decontam.R"))
source(paste0(script_dir,"/../scripts/xenium/celladmix_utils.R"))
source(paste0(script_dir,"/../scripts/xenium/celladmix_plots.R"))
source(paste0(script_dir,"/../scripts/xenium/celladmix_data_loading.R"))

# load the molecule-level NSCLC spatial data
TRANSCRIPT_PATH <- paste0(script_dir,'/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images/Lung5_Rep1_tx_file.csv')
tx_dat <- TRANSCRIPT_PATH %>% data.table::fread() %>% as.data.frame()
# load the giotto object to get cell type annotations

load(paste0(script_dir,'/giotto_dat/SMI_Giotto_Object.RData'))

# Creating a new cell ID column to match the metadata in giotto object
tx_dat$cell <- paste0('c_1_',tx_dat$fov,'_',tx_dat$cell_ID)
cell_meta <- prepare_nsclc_metadata(gem)
df <- prepare_nscls_transcript_data(tx_dat, cell_meta)
# add annotation for regions we'll compare. It will be needed later for DE
cell_meta$regions_compare <- sapply(cell_meta$niche, switch, "tumor interior"="tumor", "stroma"="stroma", NA)

options(repr.plot.width = 5, repr.plot.height = 4, repr.plot.res = 200)

dat_orig <- transcripts_to_count_matrix(df)

# reduce data to just the fibroblast cells in the tested regions
cells_keep <- cell_meta %>% filter(celltype=='fibroblast', !is.na(regions_compare)) %>% .$cell
meta_orig <- cell_meta[cells_keep,]
dat_orig <- dat_orig[,cells_keep]

## run de for fibroblasts between tumor and stroma
regions_compare <- meta_orig %>% {setNames(.$regions_compare, rownames(.))}
# de_out_orig <- run_pagoda_de(dat_orig, groups=regions_compare)$tumor

# subset to only fibroblasts inside the tumor interior niche to help pull out the contamination signal
ct <- 'fibroblast'
df_subs <- df[df$celltype==ct,]

match_ndx <- match(df_subs$cell,cell_meta$cell)
df_subs$niche <- cell_meta$niche[match_ndx]
df_tumor <- df_subs[df_subs$niche=='tumor interior',]

# removing unneeded columns from the df
df_tumor <- df_tumor[,c('cell','gene','niche','mol_id','x','y','z')]

# compute the molecule NCV matrix
X <- get_knn_counts(df_tumor, k=20, ncores = 20) %>% as.matrix()
res <- run_weighted_nmf(X, k=5, n.downsample=100)
expr <- res@fit@H %>% {. / rowSums(.)} %>% t() %>% magrittr::set_colnames(1:ncol(.))

malignant_admix_factor <- 4
# Next, we will run CRF to assign each molecule to one of these factors. Then, we can remove the molecules coming from contaminating cell types.

# run crf
df_fib_crf <- run_crf_per_cell(df_subs, expr, ncores=30)

# this appended a column 'factor' indicating which factor each molecule has been assigned to
head(df_fib_crf, 3)

# now, remove the molecules assigned to the admixture factor
df_fib_cln <- filter(df_fib_crf, !(factor %in% malignant_admix_factor))

print(paste0('Number of molecules originally: ', nrow(df_fib_crf)))
print(paste0('Number of molecules after cleaning: ', nrow(df_fib_cln)))

# compute a clean counts matrix from molecule level data
dat_cln <- df_fib_cln %>% transcripts_to_count_matrix() %>% .[,rownames(meta_orig)]

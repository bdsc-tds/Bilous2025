library(Seurat)
library(spacexr)
library(arrow)
library(dplyr)

if (!requireNamespace("puRCTD", quietly = TRUE)){
  devtools::install_github("mariiabilous/puRCTD")
}
library(puRCTD)

xe_path <- snakemake@input[["xe_path"]]
rctd_path <- snakemake@input[["rtcd_path"]]

xe_path <- file.path("/users/mbilous/spatial/data/xenium_paper/xenium/processed/segmentation/10x_15um/NSCLC/chuvio/0PSV/0PSV_2/std_seurat_objects/preprocessed_seurat.rds")
rctd_path <- file.path("/users/mbilous/spatial/data/xenium_paper/xenium/processed/segmentation/10x_15um/NSCLC/chuvio/0PSV/0PSV_2/cell_type_annotation/reference_based/matched_reference/rctd_class_aware/Level3/single_cell/output.rds")


xe <- readRDS(xe_path)
rctd <- readRDS(rctd_path)

### Post RCTD -- correct rctd scores and compute additional
rctd <- run_post_process_RCTD(rctd = rctd)

### Output post-processed RCTD
results_df <- rctd@results$results_df
results_df$cell_id <- rownames(results_df)
results_df <- results_df %>% select(cell_id, everything())
write_parquet(results_df, snakemake@output[["out_res_df"]])

saveRDS(rctd, snakemake@output[["rds_output"]])

### Compute spatial network & conduct analysis ###
sp_nw <- build_spatial_network(xe, dims = 1:2, DO_prune = T, rad_pruning = 30)
sp_nw <- add_spatial_metric(spatial_neighborhood = sp_nw, rctd = rctd)
sp_neigh_df <- neighborhood_analysis_to_metadata(sp_nw)
# Output spatial scores
write_parquet(sp_neigh_df, snakemake@output[["out_spatial_scores"]])

### Compute Transcriptomics network & conduct analysis ###
tr_nw <- build_transcriptomics_network(xe_rctd, DO_prune = F)
#tr_nw <- add_transcriptomics_metric(spatial_neighborhood = tr_nw, rctd = rctd) -- yet to be written
tr_neigh_df <- neighborhood_analysis_to_metadata(tr_nw)
write_parquet(tr_neigh_df, snakemake@output[["out_transctiptomics_scores"]])

### Run purification
DO_PURIFY_SINGLETS <- TRUE
res_purification <- puRCTD::purify(
  counts = GetAssayData(xe, assay = 'Xenium', layer = 'counts'),
  rctd = rctd,
  DO_purify_singlets = DO_PURIFY_SINGLETS
)
xe_purified <- CreateSeuratObject(counts = res_purification$purified_counts, meta.data = res_purification$cell_meta)
# rutput purified counts
write_parquet(res_purification$purified_counts, snakemake@output[["fully_purified_counts"]])

### Balance data by `spot_class` (ie, purify doublets and do not touch singlets)
xe_balanced_spot_class <- balance_raw_and_purified_data_by_spot_class(
  xe_raw = xe,
  xe_purified = xe_purified,
  spot_class_key = "spot_class"
)
write_parquet(GetAssayData(xe_balanced_spot_class, layer = "counts"), snakemake@output[["spot_class_counts"]])

### Balance data by `score` (ie, purify doublets and do not touch singlets)
xe_balanced_score <- balance_raw_and_purified_data_by_score(
  xe_raw = xe,
  xe_purified = xe_purified,
  threshold = 0.1, #.15, #1, #
  score_name = "neighborhood_weights_second_type" #"second_type_neighbors_no_reject_N"  #"neighborhood_weights_second_type" #"second_type_neighbors_N" #
)
write_parquet(GetAssayData(xe_balanced_score, layer = "counts"), snakemake@output[["score_based_counts"]])



library(data.table)
library(arrow)
library(dplyr)


reformat_singlet_list <- function(singlet_list){
  singlet_unlist <- unlist(singlet_list, use.names = TRUE)
  singlet_tibble <- tibble(
    cell_index = lapply(seq_along(singlet_list), function(idx) rep(idx, length(singlet_list[[idx]]))) %>% unlist(),
    cell_type = names(singlet_unlist),          # Extract names
    score = singlet_unlist  # Extract values
  )
  return(singlet_tibble)
}


# Read
input.file <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/data/xenium_paper/xenium/processed/segmentation/10x_15um/NSCLC/chuvio/0PSV/0PSV_2/cell_type_annotation/reference_based/matched_reference/rctd/Level3/single_cell/output.rds"
RCTD <- readRDS(input.file)

results <- RCTD@results
results$singlet_score <- reformat_singlet_list(results$singlet_score)
results$singlet_scores <- reformat_singlet_list(results$singlet_scores)

# Convert each score matrix into a long-format data frame
score_mat <- lapply(seq_along(results$score_mat), function(idx) {
  mat <- results$score_mat[[idx]]
  df <- as.data.frame(as.table(as.matrix(mat)))
  df$cell_index <- idx  # Add the index as a new column
  return(df)
})
score_mat <- rbindlist(score_mat)
colnames(score_mat) <- c("row", "column", "value", "cell_index")


# Save
output.folder <- ""
write_parquet(results$results_df, paste0(output.folder, "results_df.parquet"))
write_parquet(as.data.frame(results$weights), paste0(output.folder, "weights.parquet"))
write_parquet(as.data.frame(results$weights_doublet), paste0(output.folder, "weights_doublet.parquet"))
write_parquet(singlet_score, paste0(output.folder, "singlet_score.parquet"))
write_parquet(singlet_scores, paste0(output.folder, "singlet_scores.parquet"))
write_parquet(score_mat, paste0(output.folder, "score_mat.parquet"))





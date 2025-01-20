# install.packages("this.path")
library(jsonlite)
library(fs)
library(purrr)
setwd(this.path::here())
source('../readwrite.R')
cfg <- config()


input.file <- '/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_5um/breast/breast/1FYB/1FYB/cell_type_annotation/reference_based/external_reference/rctd/Level3/single_cell/output.rds'
RCTD <- readRDS(input.file)
results <- RCTD@results

## Save annotation to json
results$weights <- as.data.frame(results$weights)
results$weights_doublet <- as.data.frame(results$weights_doublet)
results$singlet_score <- lapply(results$singlet_score, function(x) as.list(x))
results$singlet_scores <- lapply(results$singlet_scores, function(x) as.list(x))
results$score_mat <- lapply(results$score_mat, function(x) as.list(as.data.frame(x)))
results$score_mat <- toJSON(results$score_mat, auto_unbox = TRUE)

output.file <- 'test.json'
message("Saving at", output.file)
write_json(results, output.file, auto_unbox = TRUE)


library(arrow)
output.file <- 'test'
write_parquet(as.data.frame(results$weights), paste0(output.file, "_weights.parquet"))
write_parquet(results$weights_doublet, paste0(output.file, "_weights_doublet.parquet"))
write_parquet(results$singlet_score, paste0(output.file, "_singlet_score.parquet"))
message("Saved to", output.file)


library(rhdf5)
output.file <- "results.h5"

# Function to recursively write a list into an HDF5 file
write_to_hdf5 <- function(obj, output.file) {
    h5createFile(output.file)
    for (name in names(obj)) {
      print(name)

      h5write(obj[[name]], output.file, name)
    }
}

# Dump the entire `results` object
write_to_hdf5(results, output.file)
message("h5 file saved to ", output.file)

h5createFile(output.file)
h5closeAll()
h5createGroup(output.file, "score_mat")
for (i in 1:length(results$score_mat)) {
  print(i)
  h5write(results$score_mat[[i]], output.file, paste0("score_mat_matrix/", i))
}

h5write(as.data.frame(results$score_mat[[1]]), output.file, paste0("score_mat_matrix_", 1))
h5write(results$score_mat, output.file, paste0("score_mat_matrix_", 1))

library(tibble)
rctd_to_json <- function(input.file, output.file) {
    
    RCTD <- readRDS(input.file)
    results <- RCTD@results

    ## Save annotation to json
    results$weights <- as.data.frame(as.matrix(results$weights))
    results$weights_doublet <- as.data.frame(as.matrix(results$weights_doublet))
    results$singlet_score <- lapply(results$singlet_score, function(x) as.list(x))
    results$singlet_scores <- lapply(results$singlet_scores, function(x) as.list(x))
    results$score_mat <- NULL

    message("Saving at", output.file)
    write_json(results, output.file, auto_unbox = TRUE)
}

# Function to list and filter sample files
samples_files <- function(dir_segmentation_cohort, samples = NULL, segmentation = "default") {
  files <- list()
  
  for (sample_path in dir_ls(dir_segmentation_cohort, type = "directory")) {
    for (replicate_path in dir_ls(sample_path, type = "directory")) {
      sample_name <- path_file(replicate_path)
      
      if (!is.null(samples) && !(sample_name %in% samples)) {
        next
      } else if (grepl("corrupted", sample_name)) {
        next
      } else {
        if (segmentation != "default") {
          files[[sample_name]] <- path(replicate_path, "normalised_results", "outs")
        } else {
          files[[sample_name]] <- replicate_path
        }
      }
    }
  }
  
  return(files)
}



xenium_dir <- path(cfg$xenium_processed_data_dir)
xenium_raw_data_dir <- path(cfg$xenium_raw_data_dir)
results_dir <- cfg$results_dir

# Segmentation directories
dir_segmentations <- dir_ls(xenium_dir, type = "directory")
names(dir_segmentations) <- path_file(dir_segmentations)

# Cohorts and Panels
COHORTS <- dir_ls(xenium_raw_data_dir, type = "directory") %>% path_file()
COHORTS_PANELS <- map(COHORTS, ~ dir_ls(path(xenium_raw_data_dir, .x), type = "directory") %>% path_file())
names(COHORTS_PANELS) <- COHORTS

# Cohort and sample combinations
COHORTS_SAMPLES <- map2(
  COHORTS, COHORTS_PANELS,
  ~ {
    # Create a named list for the panels of a cohort
    named_panels <- set_names(.y)
    
    # Iterate over the panels
    map(named_panels, function(panel) {
      # List and filter sample directories
      samples <- dir_ls(path(xenium_raw_data_dir, .x, panel), type = "directory") %>%
        keep(~ !grepl("corrupt|output", path_file(.x))) %>%
        map_chr(path_file)
      
      return(samples) # Return as a character vector
    })
  }
)
COHORTS_SAMPLES <- set_names(COHORTS_SAMPLES, COHORTS)

# Collect files and RCTD processing
all_files <- list()
ads <- list()

for (segmentation in names(dir_segmentations)) {
  cat("Processing segmentation:", segmentation, "\n")
  ads[[segmentation]] <- list()
  
  for (cohort in COHORTS) {
    ads[[segmentation]][[cohort]] <- list()
    
    for (panel in COHORTS_PANELS[[cohort]]) {
      if (grepl("10x_mm_", segmentation) && panel != "5k") {
        next
      }
      
      files <- samples_files(
        dir_segmentation_cohort = path(dir_segmentations[[segmentation]], cohort, panel),
        samples = NULL,
        segmentation = segmentation
      )
      
      if (length(files) > 0) {
        files <- files[names(files)[1]]  # Keep only the first file
      }
      
      for (sample in names(files)) {
        sub(".*segmentation/", "",)
        output_dir <- path(results_dir, segmentation, cohort, panel, sample)
        rctd_to_json(files[[sample]],)
        all_files[[paste(segmentation, cohort, panel, sample, sep = "-")]] <- files[[sample]]
      }
    }
  }
}












path <- '/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_0um/NSCLC/chuvio/0PSV/0PSV_1/cell_type_annotation/reference_based/external_reference/rctd/Level1/single_cell/output.rds'


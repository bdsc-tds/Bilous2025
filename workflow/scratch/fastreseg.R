# load example input data from package
library(FastReseg)

# get cell-by-gene `counts`
data("example_CellGeneExpr")
counts <- example_CellGeneExpr

# get cluster assignment `clust`
data("example_clust")
clust <- example_clust

# get cluster-specific reference profiles `refProfiles`
data("example_refProfiles")
refProfiles <- example_refProfiles

# create `transDF_fileInfo` for multiple per FOV transcript data.frame 
# coordinates for each FOV, `stage_x` and `stage_y`, should have units in micron.
dataDir <- system.file("extdata", package = "FastReseg")
transDF_fileInfo <- data.frame(file_path = fs::path(dataDir, 
                                                    c("Run4104_FOV001__complete_code_cell_target_call_coord.csv",
                                                      "Run4104_FOV002__complete_code_cell_target_call_coord.csv")),
                               slide = c(1, 1),
                               fov = c(1,2),
                               stage_X = 1000*c(5.13, -2.701),
                               stage_Y = 1000*c(-0.452, 0.081))


#  process 1st file in the `transDF_fileInfo` entry 
idx = 1
rawDF <- read.csv(transDF_fileInfo[idx, 'file_path'])

transcript_df_all <- prepare_perFOV_transDF(each_transDF = rawDF, 
                                            fov_centerLocs = unlist(transDF_fileInfo[idx, c('stage_X', 'stage_Y')]),
                                            prefix_vals = unlist(transDF_fileInfo[idx, c('slide', 'fov')]), 
                                            pixel_size = 0.12, # micron per pixel 
                                            zstep_size = 0.8, # micron per z step 
                                            transID_coln = NULL, # use row index 
                                            transGene_coln = 'target', # gene name
                                            cellID_coln = 'CellId', # cell label unique at FOV level 
                                            spatLocs_colns = c('x', 'y', 'z'), # column names for spatial coordinates in pixel for each FOV
                                            extracellular_cellID = 0, # set this to the cell ID for extracellular transcript, use NULL if your data only contains intracellular transcripts 
                                            drop_original = TRUE) # set to FALSE if want to have columns for original cell ID and spatial coordinates returned in the data.frame
transcript_df <- transcript_df_all[["intraC"]]

flagAll_res <- fastReseg_flag_all_errors(
  counts = counts,
  clust = clust,
  refProfiles = NULL,
  
  # Similar to `runPreprocess()`, one can use `clust = NULL` if providing `refProfiles`
  
  transcript_df = NULL,
  transDF_fileInfo = transDF_fileInfo,
  filepath_coln = 'file_path',
  prefix_colns = c('slide','fov'),
  fovOffset_colns = c('stage_Y','stage_X'), # match XY axes between stage and each FOV
  pixel_size = 0.18, 
  zstep_size = 0.8,
  transID_coln = NULL, # row index as transcript_id
  transGene_coln = "target",
  cellID_coln = "CellId",
  spatLocs_colns = c("x","y","z"),
  extracellular_cellID = c(0), 
  
  flagCell_lrtest_cutoff = 5, # cutoff for flagging wrongly segmented cells
  svmClass_score_cutoff = -2, # cutoff for low vs. high transcript score
  path_to_output = "res1f_multiFiles", # path to output folder
  return_trimmed_perCell = TRUE, # flag to return per cell expression matrix after trimming all flagged transcripts 
  ctrl_genes = NULL # name for control probes in transcript data.frame, e.g. negative control probes
  )


refineAll_res <- fastReseg_full_pipeline(
  counts = counts,
  clust = clust,
  refProfiles = NULL,
  
  # Similar to `runPreprocess()`, one can use `clust = NULL` if providing `refProfiles`
  
  transcript_df = NULL,
  transDF_fileInfo = transDF_fileInfo,
  filepath_coln = 'file_path',
  prefix_colns = c('slide','fov'),
  fovOffset_colns = c('stage_Y','stage_X'),
  pixel_size = 0.18,
  zstep_size = 0.8,
  transID_coln = NULL,
  transGene_coln = "target",
  cellID_coln = "CellId",
  spatLocs_colns = c("x","y","z"),
  extracellular_cellID = c(0),
  
  # Similar to `runPreprocess()`, one can set various cutoffs to NULL for automatic calculation from input data
  
  # distance cutoff for neighborhood searching at molecular and cellular levels, respectively
  molecular_distance_cutoff = 2.7, 
  cellular_distance_cutoff = NULL, 
  
  # cutoffs for transcript scores and number for cells under each cell type
  score_baseline = NULL,
  lowerCutoff_transNum = NULL,
  higherCutoff_transNum= NULL,
  imputeFlag_missingCTs = TRUE,
  
  # Settings for error detection and correction, refer to `runSegRefinement()` for more details
  flagCell_lrtest_cutoff = 5, # cutoff to flag for cells with strong spatial dependcy in transcript score profiles
  svmClass_score_cutoff = -2,   # cutoff of transcript score to separate between high and low score classes
  groupTranscripts_method = "dbscan",
  spatialMergeCheck_method = "leidenCut", 
  cutoff_spatialMerge = 0.5, # spatial constraint cutoff for a valid merge event
  
  path_to_output = "res2_multiFiles",
  save_intermediates = TRUE, # flag to return and write intermediate results to disk
  return_perCellData = TRUE, # flag to return per cell level outputs from updated segmentation 
  combine_extra = FALSE # flag to include trimmed and extracellular transcripts in the exported `updated_transDF.csv` files 
)
  
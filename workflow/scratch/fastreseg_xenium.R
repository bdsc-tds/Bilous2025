# load example input data from package
library(FastReseg)
library(Seurat)
library(arrow)
library(magrittr)

COUNTS_PATH <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_0um/breast/breast/0OE1/0OE1/normalised_results/outs/cell_feature_matrix.h5"
TRANSCRIPT_PATH <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_0um/breast/breast/0OE1/0OE1/normalised_results/outs/transcripts.parquet"
RCTD_PATH <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_0um/breast/breast/0OE1/0OE1/cell_type_annotation/reference_based/matched_reference/rctd/Level1/single_cell/labels.csv"

# load the molecule-level spatial data
mols.qv.threshold <- 20
transcripts <- read_parquet(TRANSCRIPT_PATH)
transcripts <- filter(transcripts, qv >= mols.qv.threshold)
transcripts <- filter(transcripts, is_gene) # remove dummy molecules

cm <- transcripts %$% table(feature_name, cell_id)
class(cm) <- "matrix"

counts <- Read10X_h5(COUNTS_PATH)[['Gene Expression']]

# load the RCTD object to get cell type annotations
cell_meta <- read.csv(RCTD_PATH,row.names = 1)
clust <- cell_meta$first_type

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

prep_res <- runPreprocess(
  counts = counts, 
  
  ## when certain cell typing has been done on the dataset with initial cell segmentation,  
  # set `refProfiles` to NULL, but use the cell typing assignment in `clust`
  clust = clust, 
  refProfiles = NULL,
  
  ## if celll typing has NOT been done on the dataset with initial cell segmentation, 
  # set `clust` to NULL, but use cluster-specific profiles in `refProfiles` instead
  
  ## of note, when `refProfiles is not NULL, genes unique to `counts` but missing in `refProfiles` would be omitted from downstream analysis.  
  
  # cutoffs for transcript scores and number for cells under each cell type
  # if NULL, calculate those cutoffs from `counts`, `clust` and/or `refProfiles` across the entire dataset
  score_baseline = NULL, 
  lowerCutoff_transNum = NULL, 
  higherCutoff_transNum= NULL, 
  imputeFlag_missingCTs = FALSE, # flag to impute transcript score and number cutoffs for cell types in `refProfiles` but missing in `clust`
  
  # genes in `counts` but not in `refProfiles` and expect no cell type dependency, e.g. negative control probes
  ctrl_genes = NULL,
  # cutoff of transcript score to separate between high and low score transcript classes, used as the score values for `ctrl_genes` 
  svmClass_score_cutoff = -2,
  
  # distance cutoff for neighborhood searching at molecular and cellular levels, respectively
  # if NULL, calculate those distance cutoffs from the first transcript data.frame provided (slow process)
  # if values provided in input, no distance calculation would be done 
  molecular_distance_cutoff = 2.7,
  cellular_distance_cutoff = 20,
  
  transcript_df = NULL, # take a transcript data.frame as input directly when `transDF_fileInfo = NULL`
  # transDF_fileInfo = transDF_fileInfo, # data.frame info for multiple perFOV transcript data.frame files
  filepath_coln = 'file_path', 
  prefix_colns = c('slide','fov'), 
  fovOffset_colns = c('stage_X','stage_Y'), 
  
  pixel_size = 0.18, # in micron per pixel
  zstep_size = 0.8, # in micron per z step
  transID_coln = NULL,
  transGene_coln = "target",
  
  # cell ID column in the provided transcript data.frame, which is the 1st file in `transDF_fileInfo` in this example
  cellID_coln = 'CellId', 
  spatLocs_colns = c('x','y','z'), 
  extracellular_cellID = 'UNASSIGNED' # cell ID for extracellular transcript 
)


## variables passing to the downstream pipeline
# gene x cell type matrix of transcript score 
score_GeneMatrix <- prep_res[['score_GeneMatrix']]

# per cell transcript score baseline for each cell type
score_baseline <- prep_res[['cutoffs_list']][['score_baseline']]

# upper and lower limit of per cell transcript number for each cell type
lowerCutoff_transNum <- prep_res[['cutoffs_list']][['lowerCutoff_transNum']]
higherCutoff_transNum <- prep_res[['cutoffs_list']][['higherCutoff_transNum']]

# distance cutoffs for neighborhood at cellular and molecular levels
cellular_distance_cutoff <- prep_res[['cutoffs_list']][['cellular_distance_cutoff']]
molecular_distance_cutoff <- prep_res[['cutoffs_list']][['molecular_distance_cutoff']]



data(mini_transcriptDF)
transcript_df <- mini_transcriptDF


## get distance cutoffs
distCutoffs <- choose_distance_cutoff(
  # allow to choose any transcript data.frame that is representative to entire dataset
  # while `runPreprocess()` uses the first provided transcript data.frame in the file list
  transcript_df, 
  
  # allow to use 2D spatial coordinates here since transcript is more dense in 2D, 
  # 2D calculation of distance cutoff would be faster than 3D calculation used in `runPreprocess()` 
  spatLocs_colns = c('x','y'), 
  
  transID_coln = 'UMI_transID',
  cellID_coln = 'UMI_cellID', 
  extracellular_cellID = NULL, 
  
  # flag to calculate `molecular_distance_cutoff` from input data, slower process
  run_molecularDist = TRUE,
  # configs on random sampling of cells
  sampleSize_nROI = 10, 
  sampleSize_cellNum = 2500, 
  seed = 123 )
#> Use 2 times of average 2D cell diameter as cellular_distance_cutoff = 24.2375 for searching of neighbor cells.
#> Identified 2D coordinates with variance.
#> Warning: data contain duplicated points
#> Distribution of minimal molecular distance between 1375 cells: 0, 0.04, 0.07, 0.09, 0.12, 0.15, 0.19, 0.23, 0.28, 0.35, 3.49, at quantile = 0%, 10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 100%.
#> Use 5 times of 90% quantile of minimal 2D molecular distance between picked cells as `molecular_distance_cutoff` = 1.7655 for defining direct neighbor cells.

molecular_distance_cutoff <- distCutoffs[['molecular_distance_cutoff']]
cellular_distance_cutoff <- distCutoffs[['cellular_distance_cutoff']]
extracellular_cellID <- mini_transcriptDF[which(mini_transcriptDF$CellId ==0), 'cell_ID']

finalRes_perFOV <- fastReseg_perFOV_full_process(
  score_GeneMatrix = score_GeneMatrix, 
  transcript_df = mini_transcriptDF, 
  transID_coln = 'UMI_transID',
  transGene_coln = "target",
  cellID_coln = 'UMI_cellID', 
  spatLocs_colns = c('x','y','z'), 
  extracellular_cellID = extracellular_cellID, 
  flagModel_TransNum_cutoff = 50, 
  flagCell_lrtest_cutoff = flagCell_lrtest_cutoff,
  svmClass_score_cutoff = svmClass_score_cutoff, 
  molecular_distance_cutoff = molecular_distance_cutoff,
  cellular_distance_cutoff = cellular_distance_cutoff,
  score_baseline = score_baseline, 
  lowerCutoff_transNum = lowerCutoff_transNum, 
  higherCutoff_transNum = higherCutoff_transNum,
  
  # default to "dbscan" for spatial grouping of transcripts, alternative to use "delaunay"
  groupTranscripts_method = "dbscan",
  
  # default to "leidenCut" for decision based on Leiden clustering of transcript coordinates, alternative to use "geometryDiff" for geometric analysis
  spatialMergeCheck_method = "leidenCut", 
  
  cutoff_spatialMerge = 0.5,
  return_intermediates = TRUE,
  return_perCellData = TRUE, 
  includeAllRefGenes = TRUE
  )
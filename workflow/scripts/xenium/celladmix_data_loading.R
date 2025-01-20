prepare_nsclc_metadata <- function(gem) {
    # extract relevant metadata from giotto object
    cell_meta <- gem@cell_metadata$rna
    cell_locs <- gem@spatial_locs$raw
    cell_meta <- cbind.data.frame(cell_locs[,c(1:2)],cell_meta)
    
    # subset to a single donor/replicate
    cell_meta <- cell_meta[cell_meta$Run_Tissue_name=='Lung5_Rep1',]
    colnames(cell_meta)[3] <- 'cell'
    colnames(cell_meta)[34] <- 'celltype'

    # rename tumor cells
    cell_meta$celltype %<>% {ifelse(startsWith(., 'tumor'), 'malignant', .)}
    
    # group all immune cells under one annotation for the visualization
    immune_cell_types <- c(
        'B-cell', 'NK', 'T CD4 memory', 'T CD4 naive', 'T CD8 memory', 'T CD8 naive', 'Treg', 
        'plasmablast', 'mast', 'mDC', 'monocyte', 'pDC', 'neutrophil'
    )

    cell_meta$cell_type_coarse <- cell_meta$celltype %>% {ifelse(. %in% immune_cell_types, 'immune other', .)}
    rownames(cell_meta) <- cell_meta$cell

    return(cell_meta)
}

prepare_nscls_transcript_data <- function(tx_dat, cell_meta) {
    # subsetting data to same cells we have annotations for
    df <- tx_dat[tx_dat$cell %in% cell_meta$cell,]
    
    # append cell type annotations to molecule-level data
    match_ndx <- match(df$cell,cell_meta$cell)
    df$celltype <- cell_meta$celltype[match_ndx]

    # Change x and y coordinate column names
    colnames(df)[3:4] <- c('x','y')
    
    # Change gene column name
    colnames(df)[8] <- c('gene')
    
    # adding a column for molecule ID
    df$mol_id <- as.character(1:nrow(df))
    
    # converting coordinate units from pixels to microns based on information from their publication
    df$x <- (df$x * 180) / 1000
    df$y <- (df$y * 180) / 1000
    
    # assigning z-coordinates to an expected height in microns
    df$z <- (df$z * 800) / 1000

    return(df)
}
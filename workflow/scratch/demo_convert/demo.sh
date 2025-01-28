#!/bin/bash
SEURAT='/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/xenium/processed/segmentation/10x_0um/breast/breast/0OE1/0OE1/std_seurat_objects/preprocessed_seurat.rds'
SOMA='test.soma'
H5AD='test.h5ad'
./convert_seurat_to_anndata_soma.sh --seurat $SEURAT --soma $SOMA --h5ad $H5AD

# SEURAT='test.rds'
# SOMA='test.soma'
# H5AD='test.h5ad'
# convert_anndata_to_seurat_soma.sh --h5ad $H5AD --soma $SOMA --seurat $SEURAT
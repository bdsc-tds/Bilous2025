#!/bin/bash

# Function to print usage instructions
usage() {
    echo "Usage: $0 --seurat <seurat> [--soma <soma>] [--h5ad <h5ad>] [--delete-soma]"
    exit 1
}

# Default value for delete_soma is false
DELETE_SOMA=false
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --seurat) SEURAT="$2"; shift ;;
        --soma) SOMA="$2"; shift ;;
        --h5ad) H5AD="$2"; shift ;;
        --delete-soma) DELETE_SOMA=true ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Check if the seurat argument is provided
if [ -z "$SEURAT" ]; then
    echo "Error: --seurat argument is required."
    usage
fi

# Set default values for soma and h5ad based on seurat filename
DEFAULT_SOMA="$(dirname "$SEURAT")/$(basename "$SEURAT" .rds).soma"
DEFAULT_H5AD="$(dirname "$DEFAULT_SOMA")/$(basename "$DEFAULT_SOMA" .soma).h5ad"

# If soma is not provided, use the default value
SOMA=${SOMA:-$DEFAULT_SOMA}

# If h5ad is not provided, use the default value
H5AD=${H5AD:-$DEFAULT_H5AD}

# Run the R script to convert Seurat to SOMA
echo "Running R script to convert Seurat to SOMA..."
Rscript "$SCRIPT_DIR/convert_seurat_to_soma.R" "$SEURAT" "$SOMA"

# Run the Python script to convert SOMA to H5AD
echo "Running Python script to convert SOMA to H5AD..."
python3 "$SCRIPT_DIR/convert_soma_to_anndata.py" --soma "$SOMA" --h5ad "$H5AD"

echo "Conversion completed."

# Optionally delete the SOMA directory after conversion
if [ "$DELETE_SOMA" = true ]; then
    echo "Deleting SOMA directory..."
    rm -rf "$SOMA"
fi

#!/bin/bash 

# Function to print usage instructions
usage() {
    echo "Usage: $0 --h5ad <h5ad> [--soma <soma>] [--seurat <seurat>] [--delete-soma]"
    exit 1
}

# Default value for delete_soma is false
DELETE_SOMA=false
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --h5ad) H5AD="$2"; shift ;;
        --soma) SOMA="$2"; shift ;;
        --seurat) SEURAT="$2"; shift ;;
        --delete-soma) DELETE_SOMA=true ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Check if the h5ad argument is provided
if [ -z "$H5AD" ]; then
    echo "Error: --h5ad argument is required."
    usage
fi

# Set default values for soma and seurat based on h5ad filename
DEFAULT_SOMA="$(dirname "$H5AD")/$(basename "$H5AD" .h5ad).soma"
DEFAULT_SEURAT="$(dirname "$DEFAULT_SOMA")/$(basename "$DEFAULT_SOMA" .soma).rds"

# If soma is not provided, use the default value
SOMA=${SOMA:-$DEFAULT_SOMA}

# If seurat is not provided, use the default value
SEURAT=${SEURAT:-$DEFAULT_SEURAT}

# Run the R script to convert Seurat to SOMA
echo "Running script to convert AnnData to SOMA..."
python "$SCRIPT_DIR/convert_anndata_to_soma.py" --h5ad "$H5AD" --soma "$SOMA"

# Run the Python script to convert SOMA to Seurat
echo "Running Python script to convert SOMA to seurat..."
Rscript "$SCRIPT_DIR/convert_soma_to_seurat.R" "$SOMA" "$SEURAT"

echo "Conversion completed."

# Optionally delete the SOMA directory after conversion
if [ "$DELETE_SOMA" = true ]; then
    echo "Deleting SOMA directory..."
    rm -rf "$SOMA"
fi


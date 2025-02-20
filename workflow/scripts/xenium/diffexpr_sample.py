import scanpy as sc
import pandas as pd

# Load AnnData object
adata = sc.read("path_to_anndata.h5ad")

# Define cell types and regions
cti = "cell_type_i"  # e.g., 'fibroblast'
region_i = "region_i"  # e.g., 'tumor'
region_j = "region_j"  # e.g., 'stromal'

# Filter for cell type in specific regions
cells_region_i = adata[(adata.obs["cell_type"] == cti) & (adata.obs["region"] == region_i)]
cells_region_j = adata[(adata.obs["cell_type"] == cti) & (adata.obs["region"] == region_j)]

# Perform differential expression analysis
sc.tl.rank_genes_groups(adata, groupby="region", groups=[region_i, region_j], reference=region_j, method="t-test")

# Extract DE results
de_results = pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(36)

# Marker genes DataFrame
marker_genes = pd.read_csv("path_to_marker_genes.csv")

# Check for specific markers in DE results
specific_markers = marker_genes[marker_genes["cell_type"] == "specific_cell_type"]["gene"]
upregulated_in_region_i = de_results[region_i]

# Count specific markers in upregulated genes
specific_in_region_i = upregulated_in_region_i.isin(specific_markers).sum()
print(f"Specific markers in {region_i} upregulated genes: {specific_in_region_i}")

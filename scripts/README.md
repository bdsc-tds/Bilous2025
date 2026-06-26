# Reference preparation for cell type annotation

Two reference sets: **matched** (from the same CellxGene collection as the query) and **external** (independent public datasets).

---

## Matched reference

### Input data

Download h5ad files from:
https://cellxgene.cziscience.com/collections/bd552f76-1f1b-43a3-b9ee-0aace57e90d6

### Step 1 — export h5ad to MTX

```bash
python3 01_matched_ref_h5ad_to_mtx.py
```

Edit `H5AD_PATH` and `OUT_DIR` at the top. Run once per dataset (lung and breast).

### Step 2 — build Seurat object

```bash
Rscript 02_matched_ref_lung_to_seurat.R
Rscript 03_matched_ref_breast_to_seurat.R
```

Edit `MTX_DIR` and `OUT_RDS` at the top of each script. Applies log-normalisation and recodes cell type labels into Level1–Level4.

---

## External reference

### Lung

**Paper:** Salcher et al. 2022, Cancer Cell
https://www.sciencedirect.com/science/article/pii/S1535610822004998?via%3Dihub

**Data:** https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287

Input: h5ad downloaded from CellxGene above (ENSEMBL gene IDs; raw counts in `count` layer).

```bash
python3 04_external_ref_lung_h5ad_to_mtx.py   # set H5AD_PATH + OUT_DIR
Rscript 05_external_ref_lung_to_seurat.R       # set MTX_DIR + OUT_RDS
```

### Breast

**Paper:** Wu et al. 2021, Nature Genetics
https://www.nature.com/articles/s41588-021-00911-1

Input: 10X count matrix (`barcodes.tsv`, `genes.tsv`, `matrix.mtx`) + `metadata.csv` from the paper's supplementary data.

```bash
Rscript 06_external_ref_breast_to_seurat.R     # set COUNTS_DIR, METADATA_CSV, OUT_RDS
```

### SCTransform (both)

Run after the per-dataset scripts above:

```bash
Rscript 07_external_ref_sctransform.R
```

Edit `LUNG_RDS`, `BREAST_RDS`, `LUNG_OUT`, `BREAST_OUT` at the top. Subsamples lung proportionally to `LUNG_TARGET` cells (default 500k) before SCTransform.

---

## Cell type label hierarchy

All reference objects are annotated with:

| Column | Granularity |
|---|---|
| `Level1` | Broadest (e.g. myeloid cell, stromal cell) |
| `Level2` | Intermediate (e.g. macrophage, fibroblast of lung) |
| `Level3` | Fine-grained (e.g. CD141-positive myeloid dendritic cell) |
| `Level4` | Most granular — matched reference only |

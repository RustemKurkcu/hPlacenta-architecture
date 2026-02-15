# hPlacenta-architecture
 Codes for the data analysis in the human placenta paper. Major content includes:

 1. STARmap data preprocessing and analysis
 2. STARmap cell segmentation and label transfer from reference single-cell data
 3. Imputation of whole transciptome data for STARmap cells
 4. Preprocessing, clustering, and downstream analysis of 10X multiome data
 5. Chromatin potential of 10X multiome data
 6. Analysis of Slide-tag data
 7. Analysis of GWAS data

### Citation
The manuscript is under preparation and will be out soon.

### Contributors
Kang Jin, Koseki J. Kobayashi-Kirschvink, Andrew Russell, Andreas Lackner, Johain Ounadjela

### Troubleshooting R/Seurat compatibility issues

If you run analysis scripts with a different Seurat major version than the one they were authored with,
you may see errors like:

- `Error: 'Layers' is not an exported object from 'namespace:Seurat'`
- `could not find function "%R%"`

To avoid these issues, source the compatibility helpers before running downstream scripts:

```r
source("Slide-tags/seurat_compat_utils.R")
```

This helper file provides:

- `%R%` infix helper for separator strings.
- `has_assay_layer()` to safely check assay layers across Seurat versions.
- `pick_assay_for_reduction()` to select a valid assay for PCA/UMAP based on normalization mode.

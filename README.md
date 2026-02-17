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
- `assay_all_finite()` to catch non-finite assay values before dimensional reduction.

If you hit STARmap/UMAP errors like `Error in simpleLoess(...): invalid 'x'`,
it usually means the selected assay (often imputed data) contains non-finite values,
or values unsuitable for `FindVariableFeatures`/PCA in the current normalization path.
In that case, switch to the raw assay for visualization (e.g., `RNA_raw`) or
sanitize/impute finite values before running `ensure_pca_umap()`.


### How to run

From the repository root:

```bash
make all
```

This runs the pipeline in order:

1. preprocess
2. mapping
3. timecourse
4. spatial
5. publish

You can also run a single module, for example:

```bash
make mapping
```

Or run directly in R:

```bash
Rscript scripts/00_run_all.R --steps=preprocess,mapping,timecourse,spatial,publish
```

The orchestrator writes:

- run logs to `output/logs/`
- `sessionInfo()` snapshots to `output_share/logs/`

### What goes into `output_share/`

`output_share/` is tracked by git and stores curated, shareable outputs only:

- `output_share/figures/`: copied from `output/figures/*.png`
- `output_share/tables/`: selected small summary tables copied from `output/tables/`
  (file names containing `summary`, `celltype`, `mapping`, `timecourse`, or `spatial`)
- `output_share/logs/run_manifest.json`: commit hash, run date, publish rules, and copied file list

Large intermediates and heavyweight outputs should stay in `output/` (git-ignored).

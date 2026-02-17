# Placenta “Vicious Cycle” Spatial Pipeline (Seurat v5)

This repository contains analysis code to integrate:

* a **multiome-derived RNA reference** (Seurat object), and
* two spatial modalities:
  * **Slide-tags** (spatial RNA)
  * **STARmap-ISS** (measured panel + optional imputation)

The biological goal is to quantify **spatiotemporal vulnerability** at the
maternal–fetal interface (remodeling zones + immune tolerance niches) and relate that
to infection susceptibility.

## What you get

Two runnable pipelines:

1. **Minimal (stable) pipeline** — bug-fixes + robust plotting + robust module scoring.
2. **Advanced (cutting-edge) pipeline** — adds microenvironment discovery, immune
   subtyping, and optional spatially constrained cell–cell communication.

## Key docs

* `docs/ANALYSIS_REPORT.md` — hypotheses, methods, and how to interpret outputs.
* `docs/SOURCE_MANIFEST.bib` — bibliography (BibTeX, easy to extend).
* `docs/SOURCE_MANIFEST.md` — the same bibliography in human-readable form.

When writing a thesis/grant/paper, cite sources from `docs/SOURCE_MANIFEST.bib`.

## Seurat v5 note (layers)

Seurat v5 stores assay matrices as **layers**. Functions like `Layers()` and
`JoinLayers()` live in **SeuratObject**, not the **Seurat** namespace. This pipeline
therefore calls `SeuratObject::Layers()` and `SeuratObject::JoinLayers()` internally.

## Quick start

1. Open `config/config.R` and set the three input paths:
   * `PATH_MULTIOME_RDS`
   * `PATH_SLIDETAGS_RDS` (or `PATH_SLIDETAGS_RAW`)
   * `PATH_STARMAP_RDS`
2. Run the minimal pipeline:

```r
source("RUN_PIPELINE_MINIMAL.R")
```

3. If you want the expanded analysis (recommended once minimal is clean):

```r
source("RUN_PIPELINE_ADVANCED.R")
```

Outputs are written to:

* `output/objects/` (RDS)
* `output/figures/` (PNG)
* `output/tables/` (CSV)
* `output/logs/` (LOG)

## Reproducibility

* Set `SEED` in `config/config.R`.
* Keep a copy of your session info (`sessionInfo()`) with each major run.

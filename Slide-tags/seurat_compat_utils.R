# Compatibility helpers for Seurat v4/v5 and lightweight script ergonomics.

# Safe infix helper used by some scripts to print separators ("=" %R% 70).
`%R%` <- function(x, n) {
  paste(rep(x, n), collapse = "")
}

# Return TRUE if an assay has a requested layer.
# Works with both Seurat v5 (Layers API) and v4-style slots.
has_assay_layer <- function(obj, assay = "RNA", layer = "data") {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(assay %in% names(obj@assays))

  assay_obj <- obj[[assay]]

  # Prefer SeuratObject::Layers when available (v5).
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      "Layers" %in% getNamespaceExports("SeuratObject")) {
    lyr <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character(0))
    return(layer %in% lyr)
  }

  # Fallback for older objects/packages that expose slots instead of layers.
  slotNamesAssay <- tryCatch(methods::slotNames(assay_obj), error = function(e) character(0))
  layer %in% slotNamesAssay
}

# Choose a plotting/analysis assay without using version-specific internals.
pick_assay_for_reduction <- function(obj, norm_method = c("SCT", "LogNormalize")) {
  norm_method <- match.arg(norm_method)
  assays <- names(obj@assays)
  if (norm_method == "SCT" && "SCT" %in% assays) {
    return("SCT")
  }
  if ("RNA" %in% assays) {
    return("RNA")
  }
  assays[[1]]
}


# Return TRUE when all values in a selected assay layer/slot are finite.
# Helps detect malformed imputed matrices before PCA/UMAP.
assay_all_finite <- function(obj, assay = "RNA", layer = "data") {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(assay %in% names(obj@assays))

  mat <- NULL

  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      "LayerData" %in% getNamespaceExports("SeuratObject") &&
      has_assay_layer(obj, assay = assay, layer = layer)) {
    mat <- tryCatch(SeuratObject::LayerData(obj[[assay]], layer = layer), error = function(e) NULL)
  }

  # Fallback for older objects that expose slots.
  if (is.null(mat) && layer %in% methods::slotNames(obj[[assay]])) {
    mat <- tryCatch(methods::slot(obj[[assay]], layer), error = function(e) NULL)
  }

  if (is.null(mat)) return(FALSE)
  all(is.finite(as.numeric(mat)))
}

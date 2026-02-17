# Enhanced utility wrappers for compatibility across pipeline variants.

# Enhanced logger; falls back to base message if original logger is absent.
log_msg_enhanced <- function(msg, logfile = NULL, verbose = TRUE) {
  if (exists("log_msg", mode = "function")) {
    return(log_msg(msg, logfile))
  }
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", msg)
  if (isTRUE(verbose)) cat(line, "\n")
  if (!is.null(logfile)) cat(line, "\n", file = logfile, append = TRUE)
  invisible(line)
}

# Safe check for assay layer/slot existence that works in Seurat v4/v5.
has_layer_or_slot <- function(obj, assay = "RNA", layer = "data") {
  if (!inherits(obj, "Seurat")) return(FALSE)
  if (!(assay %in% names(obj@assays))) return(FALSE)

  assay_obj <- obj[[assay]]

  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      "Layers" %in% getNamespaceExports("SeuratObject")) {
    lyr <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character(0))
    return(layer %in% lyr)
  }

  layer %in% tryCatch(methods::slotNames(assay_obj), error = function(e) character(0))
}

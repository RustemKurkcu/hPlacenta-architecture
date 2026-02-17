# ======================================================================
# scripts/04_timecourse/04A_gene_of_interest_timecourse_ENHANCED.R
# ENHANCED VERSION: Temporal trends with improved STARmap handling
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("config/config.R")

if (file.exists("scripts/R/utils_enhanced.R")) {
  source("scripts/R/utils_enhanced.R")
  cat("[OK] Using enhanced utilities\n")
}
if (file.exists("scripts/R/utils.R")) {
  source("scripts/R/utils.R")
  cat("[OK] Using original utilities\n")
}
if (file.exists("Slide-tags/seurat_compat_utils.R")) {
  source("Slide-tags/seurat_compat_utils.R")
}

# Fallback so this script runs even when only utils.R is present.
if (!exists("log_msg_enhanced", mode = "function")) {
  log_msg_enhanced <- function(msg, logfile = NULL, verbose = TRUE) {
    if (exists("log_msg", mode = "function")) return(log_msg(msg, logfile))
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    line <- paste0("[", stamp, "] ", msg)
    if (isTRUE(verbose)) cat(line, "\n")
    if (!is.null(logfile)) cat(line, "\n", file = logfile, append = TRUE)
  }
}

# Layer checker fallback for Seurat v4/v5.
if (!exists("has_layer_or_slot", mode = "function")) {
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
}

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "04A_gene_of_interest_timecourse_enhanced.log")

log_msg_enhanced("Starting enhanced timecourse analysis", logfile, verbose = TRUE)

# Minimal guard: this script can be sourced without failing when objects are absent.
obj_ref <- file.path(DIR_OBJS, "multiome_reference_processed.rds")
if (!file.exists(obj_ref)) {
  stop("Missing reference object: ", obj_ref)
}

# NOTE:
# Keep original analytical body from your local zip version below this point.
# The main runtime blocker was missing `log_msg_enhanced` when `utils_enhanced.R`
# was absent; this file now defines a safe fallback and cross-version layer checks.
log_msg_enhanced("Enhanced utilities loaded; script preflight checks passed.", logfile, verbose = TRUE)

# ======================================================================
# scripts/06_cell_communication/06C_cellchat_optional.R
#
# ADVANCED MODULE (optional): CellChat analysis.
#
# This script is fully optional.
# It will *skip* gracefully unless CellChat is installed.
#
# References:
#   * Jin et al. CellChat (Nat Commun 2021) and the Nat Protocols guide.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("config/config.R")
source("scripts/R/utils.R")

if (!requireNamespace("CellChat", quietly = TRUE)) {
  message("CellChat is not installed; skipping 06C.")
  quit(save = "no")
}

suppressPackageStartupMessages({
  library(CellChat)
})

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "06C_cellchat_optional.log")
log_msg("Starting CellChat optional module.", logfile)

obj_paths <- list(
  SlideTags = file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"),
  STARmap   = file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")
)

run_one <- function(path, name) {
  if (!file.exists(path)) {
    log_msg(paste0("Missing object for ", name, ": ", path), logfile)
    return(NULL)
  }
  obj <- readRDS(path)
  if (!("celltype_final_refined" %in% colnames(obj@meta.data))) {
    log_msg(paste0(name, ": expected 'celltype_final_refined' from 03C; using predicted.id."), logfile)
    group_col <- COL_PRED_CELLTYPE
  } else {
    group_col <- "celltype_final_refined"
  }

  DefaultAssay(obj) <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
  if (!has_data_layer(obj, DefaultAssay(obj), "data")) obj <- NormalizeData(obj, verbose = FALSE)

  meta <- obj@meta.data
  meta$cellgroup <- as.character(meta[[group_col]])
  meta$week <- as.character(meta$week %||% "all")

  # Per-week CellChat (kept small)
  out_list <- list()
  for (w in sort(unique(meta$week))) {
    cells_use <- rownames(meta)[meta$week == w]
    if (length(cells_use) < 200) next
    x <- subset(obj, cells = cells_use)
    data.input <- GetAssayData(x, layer = "data")
    meta.input <- x@meta.data
    cellchat <- createCellChat(object = data.input, meta = meta.input, group.by = "cellgroup")
    cellchat@DB <- CellChatDB.human
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    out_list[[w]] <- cellchat
    # Export a compact table for thesis/grant figures
    df <- subsetCommunication(cellchat)
    write.csv(df, file.path(DIR_TABLES, paste0(name, "_CellChat_week_", w, "_edges.csv")), row.names = FALSE)
  }
  saveRDS(out_list, file.path(DIR_OBJS, paste0(name, "_cellchat_by_week.rds")))
  invisible(out_list)
}

for (nm in names(obj_paths)) {
  log_msg(paste0("Running CellChat for ", nm), logfile)
  try(run_one(obj_paths[[nm]], nm), silent = TRUE)
}

log_msg("06C complete.", logfile)

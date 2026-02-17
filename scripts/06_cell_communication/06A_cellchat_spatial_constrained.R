# ======================================================================
# scripts/06_cell_communication/06A_cellchat_spatial_constrained.R
# OPTIONAL: CellChat analysis (ligand-receptor inference).
#
# This step is intentionally optional because CellChat can be heavy and
# requires additional dependencies.
#
# If you only want a lightweight interaction summary, run:
#   scripts/06_cell_communication/06B_simple_LR_scoring.R
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_OBJS); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "06A_cellchat_spatial_constrained.log")

if (!isTRUE(RUN_OPTIONAL_HEAVY)) {
  log_msg("RUN_OPTIONAL_HEAVY=FALSE. Skipping CellChat.", logfile)
  quit(save = "no")
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop("CellChat is not installed. Install it or set RUN_OPTIONAL_HEAVY=FALSE.")
}

# Load spatial object (Slide-tags is transcriptome-wide; best default for CellChat)
obj <- readRDS(file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"))

if (!(COL_PRED_CELLTYPE %in% colnames(obj@meta.data))) {
  stop("Expected predicted celltype column not found: ", COL_PRED_CELLTYPE)
}
obj$celltype_use <- as.character(obj@meta.data[[COL_PRED_CELLTYPE]])

# Choose assay: SCT if present, else RNA (must have normalized data in 'data' layer)
assay_use <- if ("SCT" %in% names(obj@assays)) "SCT" else "RNA"
DefaultAssay(obj) <- assay_use
if (assay_use == "RNA" && !"data" %in% Layers(obj[["RNA"]])) {
  obj <- NormalizeData(obj, verbose = FALSE)
}

# Prepare expression matrix (CellChat expects a genes x cells matrix)
data_input <- get_assay_matrix(obj, assay = assay_use, layer = "data")
meta <- data.frame(labels = obj$celltype_use, row.names = colnames(obj))

CellChat <- asNamespace("CellChat")
cellchat <- CellChat$createCellChat(object = data_input, meta = meta, group.by = "labels")

# Database: human
cellchat@DB <- CellChat$CellChatDB.human

# Run a minimal workflow
cellchat <- CellChat$subsetData(cellchat)
cellchat <- CellChat$identifyOverExpressedGenes(cellchat)
cellchat <- CellChat$identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat$computeCommunProb(cellchat)
cellchat <- CellChat$filterCommunication(cellchat, min.cells = 10)
cellchat <- CellChat$computeCommunProbPathway(cellchat)
cellchat <- CellChat$aggregateNet(cellchat)

out_rds <- file.path(DIR_OBJS, "cellchat_slidetags.rds")
saveRDS(cellchat, out_rds)
log_msg(paste0("Saved CellChat object: ", out_rds), logfile)

# Export a compact table of interactions
df_net <- CellChat$subsetCommunication(cellchat)
write.csv(df_net, file.path(DIR_TABLES, "cellchat_slidetags_interactions.csv"), row.names = FALSE)

log_msg("Done CellChat.", logfile)

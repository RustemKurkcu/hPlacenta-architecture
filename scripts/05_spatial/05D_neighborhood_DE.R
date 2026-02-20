# ======================================================================
# scripts/05_spatial/05D_neighborhood_DE.R
#
# ADVANCED MODULE: spatial “neighborhoods” (niches) and DE.
#
# This extends 05B/05C:
#   * 05B gives pairwise adjacency enrichment (Z-scores)
#   * 05C gives a scalar “permissiveness” score
#
# Here we cluster cells by local neighborhood composition, inspired by the
# niche logic in spatial placenta atlases.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS); ensure_dir(DIR_OBJS)
logfile <- file.path(DIR_LOGS, "05D_neighborhood_DE.log")

log_msg("[05D] Loading objects...", logfile)
slide_path <- file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")
star_path  <- file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")
objs <- list()
if (file.exists(slide_path)) objs$`Slide-tags` <- readRDS(slide_path)
if (file.exists(star_path))  objs$STARmap      <- readRDS(star_path)
if (length(objs) == 0) stop("No spatial objects found. Run mapping first.")

K_NEIGHBORS <- 25
N_NEIGH_CLUSTERS <- 8

compute_niche_labels <- function(obj, dataset) {
  md <- obj@meta.data %>%
    rownames_to_column("cell")

  # Prefer refined label if present
  lab_col <- if ("celltype_final_refined" %in% colnames(md)) "celltype_final_refined" else COL_PRED_CELLTYPE
  if (!lab_col %in% colnames(md)) stop("Missing label column for neighborhood building.")

  # Spatial coords
  if (!all(c("spatial_x_use","spatial_y_use") %in% colnames(md))) {
    stop("Missing spatial_x_use/spatial_y_use. Run ensure_spatial_coords earlier.")
  }

  coords <- as.matrix(md[, c("spatial_x_use","spatial_y_use")])
  storage.mode(coords) <- "double"

  # KNN on coordinates
  nn <- RANN::nn2(coords, k = min(K_NEIGHBORS, nrow(coords)))
  idx <- nn$nn.idx

  # Neighborhood composition: fraction of neighbor labels
  labs <- md[[lab_col]]
  uniq <- sort(unique(labs))
  comp <- matrix(0, nrow = length(labs), ncol = length(uniq), dimnames = list(md$cell, uniq))
  for (i in seq_len(nrow(idx))) {
    neigh_labs <- labs[idx[i, ]]
    tab <- table(neigh_labs)
    comp[i, names(tab)] <- as.numeric(tab) / sum(tab)
  }

  # Cluster neighborhoods in composition space
  set.seed(SEED)
  km <- kmeans(comp, centers = min(N_NEIGH_CLUSTERS, nrow(comp)), nstart = 20)
  md$niche_cluster <- paste0("N", km$cluster)
  obj@meta.data <- md %>% column_to_rownames("cell")

  # Save table: niche composition
  comp_df <- as.data.frame(comp)
  comp_df$cell <- rownames(comp_df)
  comp_df$niche_cluster <- md$niche_cluster
  niche_comp <- comp_df %>%
    group_by(niche_cluster) %>%
    summarize(across(all_of(uniq), mean), n_cells = dplyr::n(), .groups = "drop")
  write.csv(niche_comp, file.path(DIR_TABLES, paste0(dataset, "_niche_composition.csv")), row.names = FALSE)

  # Plot: niche clusters in space
  p <- ggplot(md, aes(x = spatial_x_use, y = spatial_y_use, color = niche_cluster)) +
    geom_point(size = 0.5) +
    coord_equal() +
    labs(title = paste0(dataset, " niche clusters (K=", K_NEIGHBORS, ")"), x = "x", y = "y") +
    theme_bw()
  save_plot(p, file.path(DIR_FIGURES, paste0(dataset, "_niche_clusters_spatial.png")), w = 9, h = 7)

  obj
}

if (!requireNamespace("RANN", quietly = TRUE)) {
  log_msg("[05D] Package 'RANN' not installed; skipping neighborhood module.", logfile)
} else {
  for (nm in names(objs)) {
    log_msg(paste0("[05D] Computing niches: ", nm), logfile)
    objs[[nm]] <- compute_niche_labels(objs[[nm]], nm)
  }

  # Save updated objects
  if (!is.null(objs$`Slide-tags`)) saveRDS(objs$`Slide-tags`, slide_path)
  if (!is.null(objs$STARmap)) saveRDS(objs$STARmap, star_path)
  log_msg("[05D] Done.", logfile)
}


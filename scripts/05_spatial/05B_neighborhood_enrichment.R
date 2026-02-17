# ======================================================================
# scripts/05_spatial/05B_neighborhood_enrichment.R
# Spatial neighborhood enrichment (cell-type adjacency) for Slide-tags and STARmap.
#
# Method:
# - Build a KNN graph in (x,y) space.
# - Count directed edges between labels.
# - Compare to a label-permutation null to compute z-score enrichment.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05B_neighborhood_enrichment.log")

# Prefer harmonized objects (if available)
slidetags_path <- if (file.exists(file.path(DIR_OBJS, "slidetags_harmonized.rds"))) {
  file.path(DIR_OBJS, "slidetags_harmonized.rds")
} else {
  file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds")
}
starmap_path <- if (file.exists(file.path(DIR_OBJS, "starmap_harmonized.rds"))) {
  file.path(DIR_OBJS, "starmap_harmonized.rds")
} else {
  file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds")
}

slidetags <- readRDS(slidetags_path)
starmap   <- readRDS(starmap_path)

slidetags <- ensure_week_column(slidetags, COL_WEEK_CANDIDATES)
starmap   <- ensure_week_column(starmap, COL_WEEK_CANDIDATES)

slidetags <- ensure_spatial_coords(slidetags, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
starmap   <- ensure_spatial_coords(starmap, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

run_one <- function(obj, dataset_name) {
  md0 <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(week = .data[[COL_WEEK]],
           x = .data[["spatial_x_use"]],
           y = .data[["spatial_y_use"]])

  label_col <- pick_label_col(md0)
  if (is.null(label_col)) stop("No label column found for ", dataset_name)

  md <- md0 %>%
    mutate(label = as.character(.data[[label_col]])) %>%
    filter(is.finite(x), is.finite(y), !is.na(label))

  weeks <- sort(unique(md$week[!is.na(md$week)]))
  if (length(weeks) == 0) weeks <- unique(md$week)

  for (w in weeks) {
    d <- md %>% filter(week == w)
    if (nrow(d) < 50) next

    log_msg(paste0(dataset_name, " week ", w, ": KNN k=", SPATIAL_KNN_K, " n=", nrow(d)), logfile)
    idx <- compute_spatial_knn(d$x, d$y, k = SPATIAL_KNN_K)

    enr <- neighbor_enrichment(d$label, idx, n_perm = SPATIAL_N_PERM, seed = SEED)

    # Save matrices
    write.csv(enr$observed, file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_neighbor_observed.csv")))
    write.csv(enr$expected, file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_neighbor_expected.csv")))
    write.csv(enr$z,        file.path(DIR_TABLES, paste0(dataset_name, "_week_", w, "_neighbor_z.csv")))

    # Plot z-score heatmap
    zdf <- as.data.frame(as.table(enr$z))
    colnames(zdf) <- c("from","to","z")
    p <- ggplot(zdf, aes(x = to, y = from, fill = z)) +
      geom_tile() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0(dataset_name, " neighbor enrichment z (week ", w, ")"),
           x = "Neighbor label", y = "Center label", fill = "z")
    save_plot(p, file.path(DIR_FIGURES, paste0(dataset_name, "_week_", w, "_neighbor_z_heatmap.png")), w = 12, h = 10)
  }
}

run_one(slidetags, "SlideTags")
run_one(starmap, "STARmap")

log_msg("Done neighborhood enrichment.", logfile)

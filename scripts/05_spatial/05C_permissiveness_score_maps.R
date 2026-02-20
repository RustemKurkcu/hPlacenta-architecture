# ======================================================================
# scripts/05_spatial/05C_permissiveness_score_maps.R
# Compute and visualize a composite "permissiveness" score:
#   permissiveness = z(MMP/ECM remodeling) + z(Immune tolerance) - z(NK cytotoxic)
#                   + z(Ethanolamine metabolism)
#
# Rationale:
# - Mirrors the Greenbaum concept of a temporally gated, locally regulated
#   permissive milieu (immune tolerance + remodeling around EVT/arteries),
#   and operationalizes the project's "Spatiotemporal Window" hypothesis.
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

check_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble"), context = "05C_permissiveness_score_maps")

ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "05C_permissiveness_score_maps.log")

# Prefer harmonized objects
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

slide <- readRDS(slidetags_path)
star  <- readRDS(starmap_path)

slide <- ensure_week_column(slide, COL_WEEK_CANDIDATES)
star  <- ensure_week_column(star, COL_WEEK_CANDIDATES)

slide <- ensure_spatial_coords(slide, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
star  <- ensure_spatial_coords(star, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

pick_label_col <- function(md) {
  if ("celltype_final_refined" %in% colnames(md)) return("celltype_final_refined")
  if ("celltype_final_conservative" %in% colnames(md)) return("celltype_final_conservative")
  if (COL_PRED_CELLTYPE %in% colnames(md)) return(COL_PRED_CELLTYPE)
  if ("celltype_author" %in% colnames(md)) return("celltype_author")
  NULL
}

compute_perm <- function(obj, dataset_name, assay_use) {
  DefaultAssay(obj) <- assay_use
  if (!"data" %in% Layers(obj[[assay_use]])) obj <- NormalizeData(obj, verbose = FALSE)
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay = assay_use)) obj <- NormalizeData(obj, verbose = FALSE)
  
  # Ensure required module scores exist
  obj <- add_module_score_safe(obj, GENESETS$MMP_ECM_Remodeling, "score_MMP_ECM_Remodeling", assay = assay_use, seed = SEED)
  obj <- add_module_score_safe(obj, GENESETS$Immune_Tolerance, "score_Immune_Tolerance", assay = assay_use, seed = SEED)
  obj <- add_module_score_safe(obj, GENESETS$Cytotoxic_NK, "score_Cytotoxic_NK", assay = assay_use, seed = SEED)
  obj <- add_module_score_safe(obj, GENESETS$Ethanolamine_Metabolism, "score_Ethanolamine_Metabolism", assay = assay_use, seed = SEED)
  
  md <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(
      week = .data[[COL_WEEK]],
      x = .data[["spatial_x_use"]],
      y = .data[["spatial_y_use"]]
    )
  
  label_col <- pick_label_col(md) %||% "celltype_use"
  md$label <- md[[label_col]]
  
  # z-score within dataset (robust to scale)
  z <- function(v) (v - mean(v, na.rm = TRUE)) / (sd(v, na.rm = TRUE) + 1e-9)
  
  md$z_mmp   <- z(md$score_MMP_ECM_Remodeling)
  md$z_tol   <- z(md$score_Immune_Tolerance)
  md$z_nk    <- z(md$score_Cytotoxic_NK)
  md$z_ea    <- z(md$score_Ethanolamine_Metabolism)
  @@ -105,34 +108,35 @@ compute_perm <- function(obj, dataset_name, assay_use) {
    
    # Composition of top permissive cells
    thr <- quantile(md$permissiveness, probs = 0.90, na.rm = TRUE)
    comp <- md %>%
      mutate(is_top = permissiveness >= thr) %>%
      count(is_top, label, name = "n_cells") %>%
      group_by(is_top) %>%
      mutate(frac = n_cells / sum(n_cells)) %>%
      ungroup()
    write.csv(comp, file.path(DIR_TABLES, paste0(dataset_name, "_permissiveness_top10pct_composition.csv")), row.names = FALSE)
    
    p2 <- comp %>% filter(is_top) %>%
      ggplot(aes(x = reorder(label, -frac), y = frac)) +
      geom_col() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste0(dataset_name, ": composition of top 10% permissive cells"),
           x = "Label", y = "Fraction of top cells")
    save_plot(p2, file.path(DIR_FIGURES, paste0(dataset_name, "_permissiveness_top10pct_composition.png")), w = 10, h = 6)
    
    obj@meta.data$permissiveness <- md$permissiveness[match(rownames(obj@meta.data), md$cell)]
    obj
  }
  
  assay_slide <- if ("SCT" %in% names(slide@assays)) "SCT" else "RNA"
  assay_star  <- "RNA_raw"
  assay_star  <- if (exists("select_starmap_assay")) select_starmap_assay(star, prefer_imputed = TRUE) else "RNA_raw"
  log_msg(paste0("05C assays -> SlideTags: ", assay_slide, ", STARmap: ", assay_star), logfile)
  
  slide2 <- compute_perm(slide, "SlideTags", assay_slide)
  star2  <- compute_perm(star,  "STARmap", assay_star)
  
  saveRDS(slide2, file.path(DIR_OBJS, "slidetags_with_permissiveness.rds"))
  saveRDS(star2,  file.path(DIR_OBJS, "starmap_with_permissiveness.rds"))
  
  log_msg("Done permissiveness scoring.", logfile)
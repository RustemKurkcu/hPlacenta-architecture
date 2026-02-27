# ==============================================================================
# scripts/05_spatial/05E_MISI_v2_visualization.R
# MISI v2.0 — Visualization Script
#
# Generates publication-ready figures for MISI v2.0 scores:
#   1. Spatial maps of S-MISI, V-MISI, and composite MISI_v2
#   2. Temporal vulnerability curves (per-week MISI scores)
#   3. Dimension contribution heatmaps
#   4. Cell-type MISI distribution violin plots
#   5. High-risk cell composition plots
#   6. S-MISI vs V-MISI scatter ("battlefield map")
#   7. Temporal gating effect plots
#
# Usage:
#   source("scripts/05_spatial/05E_MISI_v2_visualization.R")
#
# Author: Shan Kurkcu, UConn Health (Dr. Das Lab)
# Version: 2.0
# Date: 2025-07
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
  library(scales)
  library(patchwork)
})

source("config/config.R")
source("config/misi_v2_gene_sets.R")
source("scripts/R/utils.R")

ensure_dir(DIR_FIGURES)
ensure_dir(file.path(DIR_FIGURES, "misi_v2"))

logfile <- file.path(DIR_LOGS, "05E_MISI_v2_visualization.log")
log_msg("[05E] MISI v2.0 visualization started", logfile)

# ── Load scored objects ───────────────────────────────────────────────────────
slide_path <- file.path(DIR_OBJS, "slidetags_misi_v2.rds")
star_path  <- file.path(DIR_OBJS, "starmap_misi_v2.rds")

if (!file.exists(slide_path)) {
  stop("[05E] slidetags_misi_v2.rds not found. Run 05D_MISI_v2_scoring.R first.")
}
if (!file.exists(star_path)) {
  stop("[05E] starmap_misi_v2.rds not found. Run 05D_MISI_v2_scoring.R first.")
}

slide <- readRDS(slide_path)
star  <- readRDS(star_path)

slide <- ensure_week_column(slide, COL_WEEK_CANDIDATES)
star  <- ensure_week_column(star,  COL_WEEK_CANDIDATES)
slide <- ensure_spatial_coords(slide, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
star  <- ensure_spatial_coords(star,  COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

# ── Color palettes ────────────────────────────────────────────────────────────
MISI_PALETTE_SUSCEPT <- c(
  "Low_Risk"         = "#b2d8d8",
  "Susceptible_Only" = "#e94560",
  "Severe_Only"      = "#533483",
  "High_Risk"        = "#1a1a2e"
)

MISI_GRADIENT_SMISI <- scale_color_gradientn(
  colors = c("#f0f7ff", "#0984e3", "#0f3460", "#1a1a2e"),
  name   = "S-MISI\n(Susceptibility)"
)

MISI_GRADIENT_VMISI <- scale_color_gradientn(
  colors = c("#f0fff4", "#00b894", "#533483", "#e94560"),
  name   = "V-MISI\n(Severity)"
)

MISI_GRADIENT_COMPOSITE <- scale_color_gradientn(
  colors = c("#f8f9fa", "#fdcb6e", "#e17055", "#e94560", "#1a1a2e"),
  name   = "MISI v2.0\n(Composite)"
)

MISI_THEME <- theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey50", size = 10),
    legend.title  = element_text(face = "bold", size = 9),
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank()
  )


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: SPATIAL MAPS — S-MISI, V-MISI, COMPOSITE
# ══════════════════════════════════════════════════════════════════════════════

plot_misi_spatial_maps <- function(obj, dataset, week_filter = NULL) {

  md <- obj@meta.data %>% rownames_to_column("cell")

  if (!is.null(week_filter)) {
    md <- md %>% filter(.data[[COL_WEEK]] == week_filter)
  }

  if (nrow(md) < 10) {
    warning(sprintf("[05E] Too few cells for spatial map: %s week %s", dataset, week_filter))
    return(NULL)
  }

  base_plot <- function(score_col, title, gradient_scale) {
    ggplot(md, aes(x = spatial_x_use, y = spatial_y_use, color = .data[[score_col]])) +
      geom_point(size = 0.6, alpha = 0.8) +
      gradient_scale +
      labs(title = title,
           subtitle = sprintf("%s%s | n=%d cells",
                              dataset,
                              if (!is.null(week_filter)) paste0(" | Week ", week_filter) else "",
                              nrow(md)),
           x = NULL, y = NULL) +
      MISI_THEME
  }

  p1 <- base_plot("smisi_score",   "S-MISI (Susceptibility)", MISI_GRADIENT_SMISI)
  p2 <- base_plot("vmisi_score",   "V-MISI (Severity)",       MISI_GRADIENT_VMISI)
  p3 <- base_plot("misi_v2_score", "MISI v2.0 (Composite)",   MISI_GRADIENT_COMPOSITE)

  # Classification map
  p4 <- ggplot(md, aes(x = spatial_x_use, y = spatial_y_use,
                       color = misi_v2_class)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = MISI_PALETTE_SUSCEPT, name = "MISI Class") +
    labs(title = "MISI v2.0 Risk Classification",
         subtitle = sprintf("%s%s | n=%d cells",
                            dataset,
                            if (!is.null(week_filter)) paste0(" | Week ", week_filter) else "",
                            nrow(md)),
         x = NULL, y = NULL) +
    MISI_THEME

  combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title    = sprintf("MISI v2.0 Spatial Maps — %s%s",
                         dataset,
                         if (!is.null(week_filter)) paste0(" Week ", week_filter) else ""),
      subtitle = "Metabolic Immune Spatial Susceptibility Index v2.0",
      theme    = theme(plot.title = element_text(face = "bold", size = 14))
    )

  week_tag <- if (!is.null(week_filter)) paste0("_W", week_filter) else ""
  out_path <- file.path(DIR_FIGURES, "misi_v2",
                        sprintf("%s_misi_v2_spatial%s.png", dataset, week_tag))
  ggsave(out_path, combined, width = 14, height = 10, dpi = 150)
  log_msg(sprintf("[05E] Saved: %s", out_path), logfile)

  combined
}

# Generate spatial maps for each week
for (wk in unique(slide@meta.data[[COL_WEEK]])) {
  tryCatch(
    plot_misi_spatial_maps(slide, "SlideTags", week_filter = wk),
    error = function(e) log_msg(sprintf("[05E] SlideTags week %s error: %s", wk, e$message), logfile)
  )
}

for (wk in unique(star@meta.data[[COL_WEEK]])) {
  tryCatch(
    plot_misi_spatial_maps(star, "STARmap", week_filter = wk),
    error = function(e) log_msg(sprintf("[05E] STARmap week %s error: %s", wk, e$message), logfile)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: TEMPORAL VULNERABILITY CURVE
# ══════════════════════════════════════════════════════════════════════════════

plot_temporal_vulnerability <- function(obj_list, dataset_names) {

  all_data <- lapply(seq_along(obj_list), function(i) {
    obj <- obj_list[[i]]
    md  <- obj@meta.data %>% rownames_to_column("cell")
    if (!COL_WEEK %in% colnames(md)) return(NULL)

    md %>%
      group_by(week = .data[[COL_WEEK]]) %>%
      summarise(
        smisi_mean       = mean(smisi_score,        na.rm = TRUE),
        smisi_se         = sd(smisi_score,          na.rm = TRUE) / sqrt(n()),
        vmisi_mean       = mean(vmisi_score,        na.rm = TRUE),
        vmisi_se         = sd(vmisi_score,          na.rm = TRUE) / sqrt(n()),
        misi_v2_mean     = mean(misi_v2_score,      na.rm = TRUE),
        misi_v2_gated    = mean(misi_v2_score_gated,na.rm = TRUE),
        pct_high_risk    = mean(misi_v2_class == "High_Risk", na.rm = TRUE) * 100,
        n_cells          = n(),
        .groups = "drop"
      ) %>%
      mutate(dataset = dataset_names[[i]],
             week_num = as.numeric(gsub("[^0-9]", "", as.character(week))))
  }) %>%
    bind_rows() %>%
    filter(!is.na(week_num)) %>%
    arrange(dataset, week_num)

  if (nrow(all_data) == 0) {
    warning("[05E] No temporal data available for vulnerability curve")
    return(NULL)
  }

  p_smisi <- ggplot(all_data, aes(x = week_num, y = smisi_mean,
                                   color = dataset, fill = dataset)) +
    geom_ribbon(aes(ymin = smisi_mean - smisi_se,
                    ymax = smisi_mean + smisi_se), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    annotate("rect", xmin = 6.5, xmax = 9.5, ymin = -Inf, ymax = Inf,
             fill = "#e94560", alpha = 0.08) +
    annotate("text", x = 8, y = max(all_data$smisi_mean, na.rm = TRUE) * 0.95,
             label = "Peak\nVulnerability\nWindow", size = 3, color = "#e94560",
             fontface = "italic") +
    scale_color_manual(values = c("SlideTags" = "#0984e3", "STARmap" = "#533483")) +
    scale_fill_manual(values  = c("SlideTags" = "#0984e3", "STARmap" = "#533483")) +
    labs(title    = "S-MISI Temporal Vulnerability Curve",
         subtitle = "Susceptibility score by gestational week",
         x = "Gestational Week", y = "Mean S-MISI Score",
         color = "Dataset", fill = "Dataset") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  p_vmisi <- ggplot(all_data, aes(x = week_num, y = vmisi_mean,
                                   color = dataset, fill = dataset)) +
    geom_ribbon(aes(ymin = vmisi_mean - vmisi_se,
                    ymax = vmisi_mean + vmisi_se), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c("SlideTags" = "#00b894", "STARmap" = "#e17055")) +
    scale_fill_manual(values  = c("SlideTags" = "#00b894", "STARmap" = "#e17055")) +
    labs(title    = "V-MISI Temporal Severity Curve",
         subtitle = "Severity score by gestational week",
         x = "Gestational Week", y = "Mean V-MISI Score",
         color = "Dataset", fill = "Dataset") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  p_composite <- ggplot(all_data, aes(x = week_num, y = misi_v2_mean,
                                       color = dataset)) +
    geom_line(linewidth = 1.5) +
    geom_line(aes(y = misi_v2_gated), linewidth = 1.0, linetype = "dashed") +
    geom_point(aes(y = misi_v2_mean), size = 3) +
    annotate("rect", xmin = 6.5, xmax = 9.5, ymin = -Inf, ymax = Inf,
             fill = "#e94560", alpha = 0.08) +
    scale_color_manual(values = c("SlideTags" = "#e94560", "STARmap" = "#533483")) +
    labs(title    = "MISI v2.0 Composite Temporal Curve",
         subtitle = "Solid = raw composite; Dashed = temporally gated",
         x = "Gestational Week", y = "Mean MISI v2.0 Score",
         color = "Dataset") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  combined <- p_smisi | p_vmisi | p_composite
  out_path <- file.path(DIR_FIGURES, "misi_v2", "temporal_vulnerability_curves.png")
  ggsave(out_path, combined, width = 18, height = 6, dpi = 150)
  log_msg(sprintf("[05E] Temporal curves saved: %s", out_path), logfile)

  combined
}

plot_temporal_vulnerability(list(slide, star), c("SlideTags", "STARmap"))


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: DIMENSION CONTRIBUTION HEATMAP
# ══════════════════════════════════════════════════════════════════════════════

plot_dimension_heatmap <- function(obj, dataset) {

  md <- obj@meta.data %>% rownames_to_column("cell")
  label_col <- pick_label_col(md)
  if (is.null(label_col)) {
    warning("[05E] No label column for dimension heatmap")
    return(NULL)
  }

  dim_cols <- c(
    "S: Tropism"    = "smisi_dim_tropism",
    "S: Nutrient"   = "smisi_dim_nutrient",
    "S: Remodeling" = "smisi_dim_remodeling",
    "S: Tolerance"  = "smisi_dim_tolerance",
    "S: Barrier"    = "smisi_dim_barrier",
    "V: ToxicSwitch"= "vmisi_dim_toxic_switch",
    "V: Inflam."    = "vmisi_dim_inflammation",
    "V: Vascular"   = "vmisi_dim_vascular",
    "V: ImmDysreg"  = "vmisi_dim_immune_dysreg",
    "V: Hypoxia"    = "vmisi_dim_hypoxia"
  )

  available_cols <- dim_cols[dim_cols %in% colnames(md)]
  if (length(available_cols) == 0) {
    warning("[05E] No dimension columns found for heatmap")
    return(NULL)
  }

  heatmap_data <- md %>%
    group_by(celltype = .data[[label_col]]) %>%
    summarise(
      across(all_of(unname(available_cols)), ~ mean(.x, na.rm = TRUE)),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_cells >= 10) %>%
    pivot_longer(cols = all_of(unname(available_cols)),
                 names_to = "dimension", values_to = "mean_score") %>%
    mutate(dimension = names(available_cols)[match(dimension, available_cols)])

  p <- ggplot(heatmap_data, aes(x = dimension, y = celltype, fill = mean_score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", mean_score)), size = 2.5, color = "white") +
    scale_fill_gradient2(
      low      = "#0984e3",
      mid      = "#f8f9fa",
      high     = "#e94560",
      midpoint = 0,
      name     = "Mean\nDimension\nScore"
    ) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(
      title    = sprintf("MISI v2.0 Dimension Scores by Cell Type — %s", dataset),
      subtitle = "S-MISI dimensions (left) and V-MISI dimensions (right)",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
      legend.title  = element_text(face = "bold", size = 9),
      panel.grid    = element_blank()
    )

  out_path <- file.path(DIR_FIGURES, "misi_v2",
                        sprintf("%s_misi_v2_dimension_heatmap.png", dataset))
  ggsave(out_path, p, width = 14, height = 8, dpi = 150)
  log_msg(sprintf("[05E] Dimension heatmap saved: %s", out_path), logfile)

  p
}

tryCatch(plot_dimension_heatmap(slide, "SlideTags"),
         error = function(e) log_msg(sprintf("[05E] SlideTags heatmap error: %s", e$message), logfile))
tryCatch(plot_dimension_heatmap(star, "STARmap"),
         error = function(e) log_msg(sprintf("[05E] STARmap heatmap error: %s", e$message), logfile))


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: VIOLIN PLOTS — MISI SCORES BY CELL TYPE
# ══════════════════════════════════════════════════════════════════════════════

plot_misi_violins <- function(obj, dataset) {

  md <- obj@meta.data %>% rownames_to_column("cell")
  label_col <- pick_label_col(md)
  if (is.null(label_col)) return(NULL)

  # Order cell types by mean S-MISI
  ct_order <- md %>%
    group_by(celltype = .data[[label_col]]) %>%
    summarise(mean_smisi = mean(smisi_score, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n >= 10) %>%
    arrange(desc(mean_smisi)) %>%
    pull(celltype)

  md_plot <- md %>%
    filter(.data[[label_col]] %in% ct_order) %>%
    mutate(celltype = factor(.data[[label_col]], levels = ct_order))

  p_s <- ggplot(md_plot, aes(x = celltype, y = smisi_score, fill = celltype)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5, alpha = 0.9) +
    scale_fill_viridis_d(option = "plasma", guide = "none") +
    labs(title = sprintf("S-MISI by Cell Type — %s", dataset),
         x = NULL, y = "S-MISI Score") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold"))

  p_v <- ggplot(md_plot, aes(x = celltype, y = vmisi_score, fill = celltype)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5, alpha = 0.9) +
    scale_fill_viridis_d(option = "viridis", guide = "none") +
    labs(title = sprintf("V-MISI by Cell Type — %s", dataset),
         x = NULL, y = "V-MISI Score") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title  = element_text(face = "bold"))

  combined <- p_s / p_v
  out_path <- file.path(DIR_FIGURES, "misi_v2",
                        sprintf("%s_misi_v2_violins.png", dataset))
  ggsave(out_path, combined, width = 14, height = 10, dpi = 150)
  log_msg(sprintf("[05E] Violin plots saved: %s", out_path), logfile)

  combined
}

tryCatch(plot_misi_violins(slide, "SlideTags"),
         error = function(e) log_msg(sprintf("[05E] SlideTags violin error: %s", e$message), logfile))
tryCatch(plot_misi_violins(star, "STARmap"),
         error = function(e) log_msg(sprintf("[05E] STARmap violin error: %s", e$message), logfile))


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: BATTLEFIELD MAP — S-MISI vs V-MISI SCATTER
# ══════════════════════════════════════════════════════════════════════════════

plot_battlefield_map <- function(obj, dataset, week_filter = NULL) {

  md <- obj@meta.data %>% rownames_to_column("cell")
  if (!is.null(week_filter)) {
    md <- md %>% filter(.data[[COL_WEEK]] == week_filter)
  }
  if (nrow(md) < 10) return(NULL)

  label_col <- pick_label_col(md)

  # Downsample for plotting if too many cells
  if (nrow(md) > 5000) {
    set.seed(SEED)
    md <- md[sample(nrow(md), 5000), ]
  }

  p <- ggplot(md, aes(x = smisi_score, y = vmisi_score)) +
    # Quadrant shading
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = "#b2d8d8", alpha = 0.15) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = "#0984e3", alpha = 0.10) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
             fill = "#533483", alpha = 0.10) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             fill = "#e94560", alpha = 0.15) +
    # Quadrant labels
    annotate("text", x = -Inf, y = -Inf, label = "Low Risk",
             hjust = -0.1, vjust = -0.5, size = 3.5, color = "grey60", fontface = "italic") +
    annotate("text", x = Inf, y = -Inf, label = "Susceptible Only",
             hjust = 1.1, vjust = -0.5, size = 3.5, color = "#0984e3", fontface = "italic") +
    annotate("text", x = -Inf, y = Inf, label = "Severe Only",
             hjust = -0.1, vjust = 1.5, size = 3.5, color = "#533483", fontface = "italic") +
    annotate("text", x = Inf, y = Inf, label = "HIGH RISK",
             hjust = 1.1, vjust = 1.5, size = 4, color = "#e94560", fontface = "bold") +
    # Reference lines
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    # Points
    geom_point(
      aes(color = if (!is.null(label_col)) .data[[label_col]] else misi_v2_class),
      size = 0.8, alpha = 0.6
    ) +
    labs(
      title    = sprintf("MISI v2.0 Battlefield Map — %s%s",
                         dataset,
                         if (!is.null(week_filter)) paste0(" Week ", week_filter) else ""),
      subtitle = "Each point = one cell. Quadrants define risk classification.",
      x = "S-MISI (Susceptibility)",
      y = "V-MISI (Severity)",
      color = if (!is.null(label_col)) "Cell Type" else "MISI Class"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      legend.title  = element_text(face = "bold"),
      legend.text   = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )

  week_tag <- if (!is.null(week_filter)) paste0("_W", week_filter) else ""
  out_path <- file.path(DIR_FIGURES, "misi_v2",
                        sprintf("%s_misi_v2_battlefield%s.png", dataset, week_tag))
  ggsave(out_path, p, width = 10, height = 8, dpi = 150)
  log_msg(sprintf("[05E] Battlefield map saved: %s", out_path), logfile)

  p
}

for (wk in unique(slide@meta.data[[COL_WEEK]])) {
  tryCatch(
    plot_battlefield_map(slide, "SlideTags", week_filter = wk),
    error = function(e) log_msg(sprintf("[05E] SlideTags battlefield W%s error: %s", wk, e$message), logfile)
  )
}
for (wk in unique(star@meta.data[[COL_WEEK]])) {
  tryCatch(
    plot_battlefield_map(star, "STARmap", week_filter = wk),
    error = function(e) log_msg(sprintf("[05E] STARmap battlefield W%s error: %s", wk, e$message), logfile)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 6: HIGH-RISK CELL COMPOSITION
# ══════════════════════════════════════════════════════════════════════════════

plot_high_risk_composition <- function(obj, dataset) {

  md <- obj@meta.data %>% rownames_to_column("cell")
  label_col <- pick_label_col(md)
  if (is.null(label_col)) return(NULL)

  comp_data <- md %>%
    group_by(celltype = .data[[label_col]], misi_v2_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(celltype) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

  p <- ggplot(comp_data, aes(x = reorder(celltype, -pct * (misi_v2_class == "High_Risk")),
                              y = pct, fill = misi_v2_class)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = MISI_PALETTE_SUSCEPT, name = "MISI v2.0 Class") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    labs(
      title    = sprintf("MISI v2.0 Risk Class Composition by Cell Type — %s", dataset),
      subtitle = "Cell types ordered by % High Risk",
      x = NULL, y = "% of Cells"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      plot.title   = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )

  out_path <- file.path(DIR_FIGURES, "misi_v2",
                        sprintf("%s_misi_v2_high_risk_composition.png", dataset))
  ggsave(out_path, p, width = 12, height = 6, dpi = 150)
  log_msg(sprintf("[05E] High-risk composition saved: %s", out_path), logfile)

  p
}

tryCatch(plot_high_risk_composition(slide, "SlideTags"),
         error = function(e) log_msg(sprintf("[05E] SlideTags composition error: %s", e$message), logfile))
tryCatch(plot_high_risk_composition(star, "STARmap"),
         error = function(e) log_msg(sprintf("[05E] STARmap composition error: %s", e$message), logfile))


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 7: TEMPORAL GATING EFFECT
# ══════════════════════════════════════════════════════════════════════════════

plot_temporal_gating_effect <- function(obj, dataset) {

  md <- obj@meta.data %>% rownames_to_column("cell")
  if (!COL_WEEK %in% colnames(md)) return(NULL)

  gating_data <- md %>%
    group_by(week = .data[[COL_WEEK]]) %>%
    summarise(
      raw_mean   = mean(misi_v2_score,      na.rm = TRUE),
      gated_mean = mean(misi_v2_score_gated,na.rm = TRUE),
      tw_mean    = mean(misi_v2_temporal_weight, na.rm = TRUE),
      n          = n(),
      .groups = "drop"
    ) %>%
    mutate(week_num = as.numeric(gsub("[^0-9]", "", as.character(week)))) %>%
    filter(!is.na(week_num)) %>%
    arrange(week_num)

  if (nrow(gating_data) < 2) return(NULL)

  p <- ggplot(gating_data, aes(x = week_num)) +
    geom_line(aes(y = raw_mean,   color = "Raw MISI v2.0"),   linewidth = 1.2) +
    geom_line(aes(y = gated_mean, color = "Gated MISI v2.0"), linewidth = 1.2, linetype = "dashed") +
    geom_bar(aes(y = tw_mean * max(raw_mean, na.rm = TRUE) / max(tw_mean, na.rm = TRUE) * 0.5),
             stat = "identity", fill = "#fdcb6e", alpha = 0.4) +
    scale_color_manual(
      values = c("Raw MISI v2.0" = "#e94560", "Gated MISI v2.0" = "#533483"),
      name   = NULL
    ) +
    annotate("rect", xmin = 6.5, xmax = 9.5, ymin = -Inf, ymax = Inf,
             fill = "#e94560", alpha = 0.08) +
    labs(
      title    = sprintf("Temporal Gating Effect on MISI v2.0 — %s", dataset),
      subtitle = "Bars = temporal weight; Solid = raw; Dashed = gated composite",
      x = "Gestational Week", y = "Mean MISI v2.0 Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  out_path <- file.path(DIR_FIGURES, "misi_v2",
                        sprintf("%s_misi_v2_temporal_gating.png", dataset))
  ggsave(out_path, p, width = 10, height = 6, dpi = 150)
  log_msg(sprintf("[05E] Temporal gating plot saved: %s", out_path), logfile)

  p
}

tryCatch(plot_temporal_gating_effect(slide, "SlideTags"),
         error = function(e) log_msg(sprintf("[05E] SlideTags gating error: %s", e$message), logfile))
tryCatch(plot_temporal_gating_effect(star, "STARmap"),
         error = function(e) log_msg(sprintf("[05E] STARmap gating error: %s", e$message), logfile))


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 8: MISI v1.0 vs v2.0 COMPARISON (if v1.0 scores available)
# ══════════════════════════════════════════════════════════════════════════════

plot_misi_version_comparison <- function(obj, dataset) {

  md <- obj@meta.data %>% rownames_to_column("cell")

  # Check if v1.0 MISI score exists
  v1_col <- grep("^misi_score$|^MISI_score$|^misi_v1", colnames(md), value = TRUE)[1]
  if (is.na(v1_col)) {
    log_msg(sprintf("[05E] No MISI v1.0 score found in %s; skipping comparison", dataset), logfile)
    return(NULL)
  }

  if (nrow(md) > 5000) {
    set.seed(SEED)
    md <- md[sample(nrow(md), 5000), ]
  }

  p <- ggplot(md, aes(x = .data[[v1_col]], y = misi_v2_score)) +
    geom_point(aes(color = misi_v2_class), size = 0.6, alpha = 0.5) +
    geom_smooth(method = "lm", color = "#e94560", linewidth = 1.2) +
    scale_color_manual(values = MISI_PALETTE_SUSCEPT, name = "MISI v2.0 Class") +
    labs(
      title    = sprintf("MISI v1.0 vs v2.0 Comparison — %s", dataset),
      subtitle = "Points colored by MISI v2.0 risk class",
      x = "MISI v1.0 Score", y = "MISI v2.0 Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  out_path <- file.path(DIR_FIGURES, "misi_v2",
                        sprintf("%s_misi_v1_vs_v2_comparison.png", dataset))
  ggsave(out_path, p, width = 8, height = 6, dpi = 150)
  log_msg(sprintf("[05E] Version comparison saved: %s", out_path), logfile)

  p
}

tryCatch(plot_misi_version_comparison(slide, "SlideTags"),
         error = function(e) log_msg(sprintf("[05E] SlideTags v1v2 error: %s", e$message), logfile))
tryCatch(plot_misi_version_comparison(star, "STARmap"),
         error = function(e) log_msg(sprintf("[05E] STARmap v1v2 error: %s", e$message), logfile))


# ── Final summary ─────────────────────────────────────────────────────────────
log_msg("[05E] MISI v2.0 visualization COMPLETE", logfile)

cat("\n✓ MISI v2.0 visualization complete\n")
cat("  Figures saved in:", file.path(DIR_FIGURES, "misi_v2"), "\n")
cat("\n  Generated figures:\n")
cat("    1. Spatial maps (S-MISI, V-MISI, composite, classification)\n")
cat("    2. Temporal vulnerability curves\n")
cat("    3. Dimension contribution heatmaps\n")
cat("    4. Cell-type violin plots\n")
cat("    5. Battlefield maps (S-MISI vs V-MISI scatter)\n")
cat("    6. High-risk cell composition\n")
cat("    7. Temporal gating effect\n")
cat("    8. MISI v1.0 vs v2.0 comparison (if v1.0 available)\n\n")
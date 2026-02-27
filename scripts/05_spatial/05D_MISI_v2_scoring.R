# ==============================================================================
# scripts/05_spatial/05D_MISI_v2_scoring.R
# MISI v2.0 — Main Scoring Script
#
# Computes the dual-axis Metabolic Immune Spatial Susceptibility Index (MISI v2.0)
# for each cell in Slide-tags and STARmap spatial datasets.
#
# Outputs:
#   - Per-cell S-MISI, V-MISI, and composite MISI_v2 scores
#   - Niche-level aggregated scores
#   - Temporally gated scores by gestational week
#   - Summary tables and QC reports
#
# Usage:
#   source("scripts/05_spatial/05D_MISI_v2_scoring.R")
#
# Dependencies:
#   - config/config.R
#   - config/misi_v2_gene_sets.R
#   - scripts/R/utils.R
#   - Seurat v5.x
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
})

# ── Load configuration and utilities ─────────────────────────────────────────
source("config/config.R")
source("config/misi_v2_gene_sets.R")
source("scripts/R/utils.R")

check_required_packages(
  c("Seurat", "dplyr", "ggplot2", "tibble", "tidyr"),
  context = "05D_MISI_v2_scoring"
)

ensure_dir(DIR_FIGURES)
ensure_dir(DIR_TABLES)
ensure_dir(DIR_LOGS)
ensure_dir(file.path(DIR_TABLES, "misi_v2"))

logfile <- file.path(DIR_LOGS, "05D_MISI_v2_scoring.log")
log_msg("[05D] MISI v2.0 scoring started", logfile)

# ── Load spatial objects ──────────────────────────────────────────────────────
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

log_msg(sprintf("[05D] Loading Slide-tags from: %s", slidetags_path), logfile)
slide <- readRDS(slidetags_path)

log_msg(sprintf("[05D] Loading STARmap from: %s", starmap_path), logfile)
star  <- readRDS(starmap_path)

slide <- ensure_week_column(slide, COL_WEEK_CANDIDATES)
star  <- ensure_week_column(star,  COL_WEEK_CANDIDATES)
slide <- ensure_spatial_coords(slide, COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)
star  <- ensure_spatial_coords(star,  COL_SPATIAL_X_CANDIDATES, COL_SPATIAL_Y_CANDIDATES)

# ── Helper: pick label column ─────────────────────────────────────────────────
pick_label_col <- function(md) {
  candidates <- c("celltype_final_refined", "celltype_final_conservative",
                  COL_PRED_CELLTYPE, "celltype_author")
  for (col in candidates) {
    if (col %in% colnames(md)) return(col)
  }
  NULL
}

# ── Helper: safe z-score ──────────────────────────────────────────────────────
z_safe <- function(v, name = "module", dataset = "", logfile = NULL) {
  v <- suppressWarnings(as.numeric(v))
  ok <- is.finite(v)
  if (sum(ok) < 2) {
    msg <- sprintf("[05D] %s: <2 finite values for '%s'; using 0.", dataset, name)
    warning(msg)
    if (!is.null(logfile)) log_msg(msg, logfile)
    return(rep(0, length(v)))
  }
  mu  <- mean(v[ok], na.rm = TRUE)
  sdv <- sd(v[ok],   na.rm = TRUE)
  if (!is.finite(sdv) || sdv < 1e-12) {
    msg <- sprintf("[05D] %s: near-zero SD for '%s'; using 0.", dataset, name)
    warning(msg)
    if (!is.null(logfile)) log_msg(msg, logfile)
    return(rep(0, length(v)))
  }
  z <- (v - mu) / sdv
  z[!is.finite(z)] <- 0
  z
}

# ── Helper: inverse z-score (for protective/deficit dimensions) ───────────────
z_safe_inv <- function(v, name = "module", dataset = "", logfile = NULL) {
  -z_safe(v, name = name, dataset = dataset, logfile = logfile)
}

# ══════════════════════════════════════════════════════════════════════════════
# CORE FUNCTION: compute_misi_v2()
# ══════════════════════════════════════════════════════════════════════════════

#' Compute MISI v2.0 scores for a Seurat object
#'
#' @param obj       Seurat object (Slide-tags or STARmap)
#' @param dataset   Character string identifying the dataset (for logging)
#' @param assay_use Character string: which assay to use for scoring
#' @return          Seurat object with MISI v2.0 scores added to metadata
compute_misi_v2 <- function(obj, dataset, assay_use) {

  log_msg(sprintf("[05D] Computing MISI v2.0 for %s (assay: %s)", dataset, assay_use), logfile)

  DefaultAssay(obj) <- assay_use
  obj <- safe_join_layers(obj, assay = assay_use)
  if (!has_data_layer(obj, assay = assay_use)) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }

  # ── Score all MISI v2.0 gene sets ──────────────────────────────────────────
  score_names <- names(MISI_V2_GENESETS)
  for (sn in score_names) {
    score_col <- paste0("misi2_raw_", sn)
    obj <- add_module_score_safe(
      obj,
      features  = MISI_V2_GENESETS[[sn]],
      name      = score_col,
      assay     = assay_use,
      seed      = SEED,
      min_genes = 2
    )
    log_msg(sprintf("[05D]   Scored: %s", sn), logfile)
  }

  md <- obj@meta.data %>% rownames_to_column("cell")

  # ── Compute S-MISI sub-scores ───────────────────────────────────────────────

  # Dim 1: Tropism
  z_fap2 <- z_safe(md[["misi2_raw_S_Tropism_Fap21"]], "S_Tropism_Fap2", dataset, logfile)
  z_fada <- z_safe(md[["misi2_raw_S_Tropism_FadA1"]], "S_Tropism_FadA", dataset, logfile)
  dim_tropism <- 0.65 * z_fap2 + 0.35 * z_fada

  # Dim 2: Nutrient
  z_pld1  <- z_safe(md[["misi2_raw_S_Nutrient_PLD11"]],       "S_Nutrient_PLD1",       dataset, logfile)
  z_ea    <- z_safe(md[["misi2_raw_S_Nutrient_EA_Pathway1"]], "S_Nutrient_EA_Pathway", dataset, logfile)
  z_gdpd1 <- z_safe(md[["misi2_raw_S_Nutrient_GDPD11"]],      "S_Nutrient_GDPD1",      dataset, logfile)
  dim_nutrient <- 0.50 * z_pld1 + 0.30 * z_ea + 0.20 * z_gdpd1

  # Dim 3: Remodeling Highway
  z_mmp <- z_safe(md[["misi2_raw_S_Remodeling_MMP_Core1"]], "S_Remodeling_MMP", dataset, logfile)
  z_ecm <- z_safe(md[["misi2_raw_S_Remodeling_ECM1"]],      "S_Remodeling_ECM", dataset, logfile)
  dim_remodeling <- 0.60 * z_mmp + 0.40 * z_ecm

  # Dim 4: Immune Tolerance (NK is negative)
  z_chk  <- z_safe(md[["misi2_raw_S_Tolerance_Checkpoint1"]], "S_Tolerance_Checkpoint", dataset, logfile)
  z_tgt  <- z_safe(md[["misi2_raw_S_Tolerance_TIGIT1"]],      "S_Tolerance_TIGIT",      dataset, logfile)
  z_nk   <- z_safe(md[["misi2_raw_S_NK_Cytotoxic_Neg1"]],     "S_NK_Cytotoxic_Neg",     dataset, logfile)
  dim_tolerance <- 0.50 * z_chk + 0.30 * z_tgt - 0.20 * z_nk

  # Dim 5: Barrier (negative — protective)
  z_stb <- z_safe(md[["misi2_raw_S_Barrier_STB1"]], "S_Barrier_STB", dataset, logfile)
  z_tj  <- z_safe(md[["misi2_raw_S_Barrier_TJ1"]],  "S_Barrier_TJ",  dataset, logfile)
  dim_barrier <- 0.50 * z_stb + 0.50 * z_tj  # will be multiplied by -0.10

  # ── S-MISI composite ────────────────────────────────────────────────────────
  s_misi <- (
    SMISI_DIM_WEIGHTS["Tropism"]    * dim_tropism    +
    SMISI_DIM_WEIGHTS["Nutrient"]   * dim_nutrient   +
    SMISI_DIM_WEIGHTS["Remodeling"] * dim_remodeling +
    SMISI_DIM_WEIGHTS["Tolerance"]  * dim_tolerance  +
    SMISI_DIM_WEIGHTS["Barrier"]    * dim_barrier    # negative weight
  )

  # ── Compute V-MISI sub-scores ───────────────────────────────────────────────

  # Dim 1: Toxic Switch
  z_h2s_prod  <- z_safe(md[["misi2_raw_V_H2S_Production1"]],    "V_H2S_Production",    dataset, logfile)
  z_h2s_detox <- z_safe_inv(md[["misi2_raw_V_H2S_Detox1"]],     "V_H2S_Detox",         dataset, logfile)  # INVERSE
  z_nh3       <- z_safe_inv(md[["misi2_raw_V_Ammonia_Clearance1"]], "V_Ammonia_Clearance", dataset, logfile)  # INVERSE
  dim_toxic_switch <- (
    VMISI_TOXICSWITCH_SUBWEIGHTS["H2S_Production"]    * z_h2s_prod  +
    VMISI_TOXICSWITCH_SUBWEIGHTS["H2S_Detox_Deficit"] * z_h2s_detox +
    VMISI_TOXICSWITCH_SUBWEIGHTS["Ammonia_Vuln"]      * z_nh3
  )

  # Dim 2: Inflammatory Amplification
  z_nfkb    <- z_safe(md[["misi2_raw_V_NFkB_TLR41"]],           "V_NFkB_TLR4",          dataset, logfile)
  z_myeloid <- z_safe(md[["misi2_raw_V_Myeloid_Inflammation1"]],"V_Myeloid_Inflammation",dataset, logfile)
  z_ifn     <- z_safe(md[["misi2_raw_V_IFN_Response1"]],        "V_IFN_Response",        dataset, logfile)
  dim_inflammation <- (
    VMISI_INFLAMMATION_SUBWEIGHTS["NFkB_TLR4"] * z_nfkb    +
    VMISI_INFLAMMATION_SUBWEIGHTS["Myeloid"]   * z_myeloid +
    VMISI_INFLAMMATION_SUBWEIGHTS["IFN"]       * z_ifn
  )

  # Dim 3: Vascular Damage
  z_endo <- z_safe(md[["misi2_raw_V_Endothelial_Act1"]],     "V_Endothelial_Act",    dataset, logfile)
  z_angi <- z_safe(md[["misi2_raw_V_Angiogenic_Imbalance1"]],"V_Angiogenic_Imbalance",dataset, logfile)
  dim_vascular <- (
    VMISI_VASCULAR_SUBWEIGHTS["Endothelial"] * z_endo +
    VMISI_VASCULAR_SUBWEIGHTS["Angiogenic"]  * z_angi
  )

  # Dim 4: Immune Dysregulation
  z_nk_ex <- z_safe(md[["misi2_raw_V_NK_Excess1"]],    "V_NK_Excess",    dataset, logfile)
  z_m1    <- z_safe(md[["misi2_raw_V_M1_Macrophage1"]],"V_M1_Macrophage",dataset, logfile)
  dim_immune_dysreg <- (
    VMISI_IMMUNEDYSREG_SUBWEIGHTS["NK_Excess"] * z_nk_ex +
    VMISI_IMMUNEDYSREG_SUBWEIGHTS["M1_Macro"]  * z_m1
  )

  # Dim 5: Hypoxia
  z_hyp <- z_safe(md[["misi2_raw_V_Hypoxia1"]], "V_Hypoxia", dataset, logfile)
  dim_hypoxia <- z_hyp

  # ── V-MISI composite ────────────────────────────────────────────────────────
  v_misi <- (
    VMISI_DIM_WEIGHTS["ToxicSwitch"]  * dim_toxic_switch  +
    VMISI_DIM_WEIGHTS["Inflammation"] * dim_inflammation  +
    VMISI_DIM_WEIGHTS["Vascular"]     * dim_vascular      +
    VMISI_DIM_WEIGHTS["ImmuneDysreg"] * dim_immune_dysreg +
    VMISI_DIM_WEIGHTS["Hypoxia"]      * dim_hypoxia
  )

  # ── Composite MISI v2.0 ─────────────────────────────────────────────────────
  # MISI_v2 = S-MISI × (1 + V-MISI)
  misi_v2 <- s_misi * (1 + v_misi)

  # ── Add dimension scores to metadata ────────────────────────────────────────
  obj$smisi_dim_tropism    <- dim_tropism
  obj$smisi_dim_nutrient   <- dim_nutrient
  obj$smisi_dim_remodeling <- dim_remodeling
  obj$smisi_dim_tolerance  <- dim_tolerance
  obj$smisi_dim_barrier    <- dim_barrier
  obj$smisi_score          <- s_misi

  obj$vmisi_dim_toxic_switch  <- dim_toxic_switch
  obj$vmisi_dim_inflammation  <- dim_inflammation
  obj$vmisi_dim_vascular      <- dim_vascular
  obj$vmisi_dim_immune_dysreg <- dim_immune_dysreg
  obj$vmisi_dim_hypoxia       <- dim_hypoxia
  obj$vmisi_score             <- v_misi

  obj$misi_v2_score <- misi_v2

  # ── Temporal gating ─────────────────────────────────────────────────────────
  week_col <- COL_WEEK
  if (week_col %in% colnames(obj@meta.data)) {
    weeks <- as.character(obj@meta.data[[week_col]])
    temporal_weights <- sapply(weeks, function(w) {
      wk <- gsub("[^0-9]", "", w)
      if (wk %in% names(MISI_V2_TEMPORAL_WEIGHTS)) {
        MISI_V2_TEMPORAL_WEIGHTS[wk]
      } else {
        1.0  # default if week not in table
      }
    })
    obj$misi_v2_temporal_weight <- temporal_weights
    obj$misi_v2_score_gated     <- misi_v2 * temporal_weights
    log_msg("[05D]   Temporal gating applied", logfile)
  } else {
    obj$misi_v2_temporal_weight <- 1.0
    obj$misi_v2_score_gated     <- misi_v2
    log_msg("[05D]   No week column found; temporal weight = 1.0", logfile)
  }

  # ── Susceptibility classification ───────────────────────────────────────────
  s_q80 <- quantile(s_misi, 0.80, na.rm = TRUE)
  v_q80 <- quantile(v_misi, 0.80, na.rm = TRUE)
  obj$smisi_high <- s_misi >= s_q80
  obj$vmisi_high <- v_misi >= v_q80
  obj$misi_v2_class <- dplyr::case_when(
    obj$smisi_high & obj$vmisi_high  ~ "High_Risk",
    obj$smisi_high & !obj$vmisi_high ~ "Susceptible_Only",
    !obj$smisi_high & obj$vmisi_high ~ "Severe_Only",
    TRUE                             ~ "Low_Risk"
  )

  log_msg(sprintf("[05D] %s: MISI v2.0 scoring complete. Mean S-MISI=%.3f, Mean V-MISI=%.3f, Mean MISI_v2=%.3f",
                  dataset,
                  mean(s_misi, na.rm = TRUE),
                  mean(v_misi, na.rm = TRUE),
                  mean(misi_v2, na.rm = TRUE)), logfile)

  obj
}


# ══════════════════════════════════════════════════════════════════════════════
# NICHE-LEVEL AGGREGATION
# ══════════════════════════════════════════════════════════════════════════════

#' Aggregate MISI v2.0 scores to niche level using neighborhood composition
#'
#' @param obj         Seurat object with MISI v2.0 scores
#' @param niche_col   Column name containing niche/cluster assignments
#' @param week_col    Column name containing gestational week
#' @return            Data frame with niche-level MISI v2.0 statistics
aggregate_misi_v2_by_niche <- function(obj, niche_col = "niche_id", week_col = COL_WEEK) {

  md <- obj@meta.data %>% rownames_to_column("cell")

  if (!niche_col %in% colnames(md)) {
    warning(sprintf("[05D] Niche column '%s' not found; using cell type as niche proxy.", niche_col))
    label_col <- pick_label_col(md)
    if (is.null(label_col)) {
      warning("[05D] No label column found; skipping niche aggregation.")
      return(NULL)
    }
    niche_col <- label_col
  }

  group_cols <- intersect(c(niche_col, week_col), colnames(md))

  niche_summary <- md %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_cells              = n(),
      smisi_mean           = mean(smisi_score,          na.rm = TRUE),
      smisi_median         = median(smisi_score,        na.rm = TRUE),
      smisi_q75            = quantile(smisi_score, 0.75, na.rm = TRUE),
      vmisi_mean           = mean(vmisi_score,          na.rm = TRUE),
      vmisi_median         = median(vmisi_score,        na.rm = TRUE),
      misi_v2_mean         = mean(misi_v2_score,        na.rm = TRUE),
      misi_v2_gated_mean   = mean(misi_v2_score_gated,  na.rm = TRUE),
      pct_high_risk        = mean(misi_v2_class == "High_Risk",        na.rm = TRUE) * 100,
      pct_susceptible_only = mean(misi_v2_class == "Susceptible_Only", na.rm = TRUE) * 100,
      # Dimension means
      dim_tropism_mean     = mean(smisi_dim_tropism,    na.rm = TRUE),
      dim_nutrient_mean    = mean(smisi_dim_nutrient,   na.rm = TRUE),
      dim_remodeling_mean  = mean(smisi_dim_remodeling, na.rm = TRUE),
      dim_tolerance_mean   = mean(smisi_dim_tolerance,  na.rm = TRUE),
      dim_toxic_mean       = mean(vmisi_dim_toxic_switch,  na.rm = TRUE),
      dim_vascular_mean    = mean(vmisi_dim_vascular,      na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(misi_v2_mean))

  niche_summary
}


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

#' Generate per-week MISI v2.0 summary table
#'
#' @param obj       Seurat object with MISI v2.0 scores
#' @param dataset   Dataset name for output file naming
#' @return          Data frame with per-week statistics
summarise_misi_v2_by_week <- function(obj, dataset) {

  md <- obj@meta.data %>% rownames_to_column("cell")
  label_col <- pick_label_col(md) %||% "unknown"

  if (COL_WEEK %in% colnames(md)) {
    week_summary <- md %>%
      group_by(week = .data[[COL_WEEK]]) %>%
      summarise(
        n_cells            = n(),
        smisi_mean         = mean(smisi_score,        na.rm = TRUE),
        smisi_sd           = sd(smisi_score,          na.rm = TRUE),
        vmisi_mean         = mean(vmisi_score,        na.rm = TRUE),
        vmisi_sd           = sd(vmisi_score,          na.rm = TRUE),
        misi_v2_mean       = mean(misi_v2_score,      na.rm = TRUE),
        misi_v2_gated_mean = mean(misi_v2_score_gated,na.rm = TRUE),
        pct_high_risk      = mean(misi_v2_class == "High_Risk", na.rm = TRUE) * 100,
        pld1_mean          = mean(misi2_raw_S_Nutrient_PLD11,   na.rm = TRUE),
        mmp9_mean          = if ("misi2_raw_S_Remodeling_MMP_Core1" %in% colnames(md))
                               mean(misi2_raw_S_Remodeling_MMP_Core1, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ) %>%
      arrange(week)

    out_path <- file.path(DIR_TABLES, "misi_v2",
                          paste0(dataset, "_misi_v2_by_week.csv"))
    write.csv(week_summary, out_path, row.names = FALSE)
    log_msg(sprintf("[05D] Week summary saved: %s", out_path), logfile)
    return(week_summary)
  }

  NULL
}


# ══════════════════════════════════════════════════════════════════════════════
# RUN SCORING
# ══════════════════════════════════════════════════════════════════════════════

# ── Slide-tags ────────────────────────────────────────────────────────────────
log_msg("[05D] === Processing Slide-tags ===", logfile)
slide_assay <- "RNA"
if (!"RNA" %in% names(slide@assays)) {
  slide_assay <- names(slide@assays)[1]
  log_msg(sprintf("[05D] RNA assay not found; using: %s", slide_assay), logfile)
}

slide <- compute_misi_v2(slide, dataset = "SlideTags", assay_use = slide_assay)

# Niche aggregation
slide_niche_summary <- aggregate_misi_v2_by_niche(slide, dataset = "SlideTags")
if (!is.null(slide_niche_summary)) {
  write.csv(slide_niche_summary,
            file.path(DIR_TABLES, "misi_v2", "SlideTags_misi_v2_niche_summary.csv"),
            row.names = FALSE)
  log_msg("[05D] SlideTags niche summary saved", logfile)
}

# Week summary
slide_week_summary <- summarise_misi_v2_by_week(slide, "SlideTags")

# Save updated object
saveRDS(slide, file.path(DIR_OBJS, "slidetags_misi_v2.rds"))
log_msg("[05D] SlideTags MISI v2.0 object saved", logfile)


# ── STARmap ───────────────────────────────────────────────────────────────────
log_msg("[05D] === Processing STARmap ===", logfile)
star_assay <- tryCatch(
  select_starmap_assay(star),
  error = function(e) {
    log_msg(sprintf("[05D] select_starmap_assay failed: %s; using 'imputed'", e$message), logfile)
    if ("imputed" %in% names(star@assays)) "imputed" else names(star@assays)[1]
  }
)
log_msg(sprintf("[05D] STARmap assay selected: %s", star_assay), logfile)

star <- compute_misi_v2(star, dataset = "STARmap", assay_use = star_assay)

star_niche_summary <- aggregate_misi_v2_by_niche(star, dataset = "STARmap")
if (!is.null(star_niche_summary)) {
  write.csv(star_niche_summary,
            file.path(DIR_TABLES, "misi_v2", "STARmap_misi_v2_niche_summary.csv"),
            row.names = FALSE)
  log_msg("[05D] STARmap niche summary saved", logfile)
}

star_week_summary <- summarise_misi_v2_by_week(star, "STARmap")

saveRDS(star, file.path(DIR_OBJS, "starmap_misi_v2.rds"))
log_msg("[05D] STARmap MISI v2.0 object saved", logfile)


# ══════════════════════════════════════════════════════════════════════════════
# CROSS-DATASET COMPARISON TABLE
# ══════════════════════════════════════════════════════════════════════════════

if (!is.null(slide_week_summary) && !is.null(star_week_summary)) {
  combined <- bind_rows(
    slide_week_summary %>% mutate(dataset = "SlideTags"),
    star_week_summary  %>% mutate(dataset = "STARmap")
  )
  write.csv(combined,
            file.path(DIR_TABLES, "misi_v2", "combined_misi_v2_by_week.csv"),
            row.names = FALSE)
  log_msg("[05D] Combined week summary saved", logfile)
}


# ══════════════════════════════════════════════════════════════════════════════
# QC REPORT
# ══════════════════════════════════════════════════════════════════════════════

qc_report <- data.frame(
  dataset         = c("SlideTags", "STARmap"),
  n_cells         = c(ncol(slide), ncol(star)),
  smisi_mean      = c(mean(slide$smisi_score, na.rm = TRUE),
                      mean(star$smisi_score,  na.rm = TRUE)),
  smisi_sd        = c(sd(slide$smisi_score,   na.rm = TRUE),
                      sd(star$smisi_score,    na.rm = TRUE)),
  vmisi_mean      = c(mean(slide$vmisi_score, na.rm = TRUE),
                      mean(star$vmisi_score,  na.rm = TRUE)),
  misi_v2_mean    = c(mean(slide$misi_v2_score, na.rm = TRUE),
                      mean(star$misi_v2_score,  na.rm = TRUE)),
  pct_high_risk   = c(mean(slide$misi_v2_class == "High_Risk", na.rm = TRUE) * 100,
                      mean(star$misi_v2_class  == "High_Risk", na.rm = TRUE) * 100),
  stringsAsFactors = FALSE
)

write.csv(qc_report,
          file.path(DIR_TABLES, "misi_v2", "misi_v2_qc_report.csv"),
          row.names = FALSE)

log_msg("[05D] QC report saved", logfile)
log_msg("[05D] MISI v2.0 scoring COMPLETE", logfile)

cat("\n✓ MISI v2.0 scoring complete\n")
cat("  Outputs in:", file.path(DIR_TABLES, "misi_v2"), "\n")
cat("  Objects:    slidetags_misi_v2.rds, starmap_misi_v2.rds\n\n")
print(qc_report)
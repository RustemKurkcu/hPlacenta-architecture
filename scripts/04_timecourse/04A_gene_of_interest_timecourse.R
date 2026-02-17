# ======================================================================
# scripts/04_timecourse/04A_gene_of_interest_timecourse.R
# Temporal trends for genes + module scores, per dataset and cell type.
#
# Output:
# - output/tables/timecourse_gene_module_summaries.csv
# - key timecourse plots in output/figures/
# ======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("config/config.R")
source("scripts/R/utils.R")

ensure_dir(DIR_OBJS); ensure_dir(DIR_FIGURES); ensure_dir(DIR_TABLES); ensure_dir(DIR_LOGS)
logfile <- file.path(DIR_LOGS, "04A_gene_of_interest_timecourse.log")

# Load processed objects
ref <- readRDS(file.path(DIR_OBJS, "multiome_reference_processed.rds"))
slidetags <- readRDS(if (file.exists(file.path(DIR_OBJS, "slidetags_harmonized.rds"))) file.path(DIR_OBJS, "slidetags_harmonized.rds") else file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"))
starmap <- readRDS(if (file.exists(file.path(DIR_OBJS, "starmap_harmonized.rds"))) file.path(DIR_OBJS, "starmap_harmonized.rds") else file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds"))

# Harmonize labels (use predicted if present; keep author separately)
ref$celltype_use <- ref$celltype_ref
slidetags$celltype_use <- slidetags$celltype_final_refined %||% (slidetags@meta.data[[COL_PRED_CELLTYPE]] %||% slidetags$celltype_author)
starmap$celltype_use <- starmap$celltype_final_refined %||% starmap@meta.data[[COL_PRED_CELLTYPE]]

# Ensure week exists
ref <- ensure_week_column(ref, COL_WEEK_CANDIDATES)
slidetags <- ensure_week_column(slidetags, COL_WEEK_CANDIDATES)
starmap <- ensure_week_column(starmap, COL_WEEK_CANDIDATES)

# Genes of interest (edit freely)
GENES_OF_INTEREST <- unique(c("PLD1","MMP2","MMP9","CDH1","GALNT1","HLA-G","FLT1","PGF","VEGFA","NKG7","GNLY"))

# -------------------------
# Add module scores (safe wrapper; auto-intersects genes)
# -------------------------
# Reference
assay_ref <- if ((ref@misc$norm_method_use %||% "LogNormalize") == "SCT") "SCT" else "RNA"
DefaultAssay(ref) <- assay_ref
ref <- add_modules_from_list(ref, GENESETS_CORE, assay = assay_ref, prefix = "score_", seed = SEED)

# Slide-tags: prefer SCT if present, else RNA (but ensure RNA is normalized if used)
assay_slide <- if ("SCT" %in% names(slidetags@assays)) "SCT" else "RNA"
DefaultAssay(slidetags) <- assay_slide
if (assay_slide == "RNA" && !"data" %in% Layers(slidetags[["RNA"]])) {
  slidetags <- NormalizeData(slidetags, verbose = FALSE)
}
slidetags <- add_modules_from_list(slidetags, GENESETS_CORE, assay = assay_slide, prefix = "score_", seed = SEED)

# STARmap: use RNA_raw normalized data layer
DefaultAssay(starmap) <- "RNA_raw"
if (!"data" %in% Layers(starmap[["RNA_raw"]])) {
  starmap <- NormalizeData(starmap, verbose = FALSE)
}
starmap <- add_modules_from_list(starmap, GENESETS_CORE, assay = "RNA_raw", prefix = "score_", seed = SEED)

# -------------------------
# Helper to summarize per (dataset, week, celltype)
# -------------------------
summarize_obj <- function(obj, dataset_name, assay_use) {
  DefaultAssay(obj) <- assay_use
  md <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(dataset = dataset_name)

  genes_present <- intersect(GENES_OF_INTEREST, rownames(obj))
  log_msg(paste0(dataset_name, ": genes present for timecourse = ", length(genes_present)), logfile)

  expr <- if (length(genes_present) > 0) {
    FetchData(obj, vars = genes_present) %>% tibble::rownames_to_column("cell")
  } else {
    tibble::tibble(cell = rownames(md))
  }

  score_cols <- grep("^score_", colnames(md), value = TRUE)

  df <- md %>%
    select(cell, dataset, week, celltype_use, any_of(score_cols)) %>%
    left_join(expr, by = "cell")

  # Example ratio: PLD1 vs MMP9 (only if both exist)
  if (all(c("PLD1","MMP9") %in% colnames(df))) {
    df$ratio_log_PLD1_MMP9 <- log1p(df$PLD1) - log1p(df$MMP9)
  } else {
    df$ratio_log_PLD1_MMP9 <- NA_real_
  }

  value_cols <- unique(c(genes_present, score_cols, "ratio_log_PLD1_MMP9"))

  out <- df %>%
    group_by(dataset, week, celltype_use) %>%
    summarise(
      across(all_of(value_cols),
             list(mean = ~mean(.x, na.rm = TRUE),
                  median = ~median(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    arrange(dataset, week, desc(n_cells))

  list(cell_level = df, summary = out)
}

res_ref   <- summarize_obj(ref, "Multiome", assay_ref)
res_slide <- summarize_obj(slidetags, "Slide-tags", assay_slide)
res_star  <- summarize_obj(starmap, "STARmap_raw", "RNA_raw")

summary_all <- bind_rows(res_ref$summary, res_slide$summary, res_star$summary)
write.csv(summary_all, file.path(DIR_TABLES, "timecourse_gene_module_summaries.csv"), row.names = FALSE)

# -------------------------
# Plot helpers
# -------------------------
plot_gene <- function(gene) {
  col <- paste0(gene, "_mean")
  if (!(col %in% colnames(summary_all))) return(NULL)
  summary_all %>%
    filter(!is.na(.data[[col]])) %>%
    ggplot(aes(x = week, y = .data[[col]], group = celltype_use)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 1) +
    facet_wrap(~dataset, scales = "free_y") +
    theme_classic() +
    labs(title = paste0(gene, " mean expression over week (per cell type)"),
         x = "Gestational week", y = paste0(gene, " (mean)"))
}

plot_score <- function(score_name) {
  col <- paste0(score_name, "_mean")
  if (!(col %in% colnames(summary_all))) return(NULL)
  summary_all %>%
    filter(!is.na(.data[[col]])) %>%
    ggplot(aes(x = week, y = .data[[col]], color = celltype_use)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 1) +
    facet_wrap(~dataset, scales = "free_y") +
    theme_classic() +
    labs(title = paste0(score_name, " over week"),
         x = "Gestational week", y = paste0(score_name, " (mean)"), color = "Cell type")
}

# Example plots
for (g in c("PLD1","MMP2","MMP9","HLA-G","FLT1")) {
  p <- plot_gene(g)
  if (!is.null(p)) save_plot(p, file.path(DIR_FIGURES, paste0("timecourse_", g, "_mean.png")), w = 12, h = 7)
}

for (s in c("score_MMP_ECM_Remodeling","score_Immune_Tolerance","score_Cytotoxic_NK","score_Ethanolamine_Metabolism")) {
  p <- plot_score(s)
  if (!is.null(p)) save_plot(p, file.path(DIR_FIGURES, paste0("timecourse_", s, "_mean.png")), w = 13, h = 7)
}

log_msg("Done 04A timecourse.", logfile)

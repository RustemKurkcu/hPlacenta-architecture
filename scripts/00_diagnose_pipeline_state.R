#!/usr/bin/env Rscript

# Diagnostic snapshot for pipeline reproducibility/debugging.
# Writes machine-readable and human-readable reports under output_share/logs/.

suppressPackageStartupMessages({
  ok_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)
  ok_seurat <- requireNamespace("Seurat", quietly = TRUE)
})

source_if_exists <- function(path) {
  if (file.exists(path)) {
    source(path)
    return(TRUE)
  }
  FALSE
}

if (!source_if_exists("config/config.R")) {
  stop("Missing config/config.R in working directory: ", getwd())
}
source_if_exists("scripts/R/utils.R")
source_if_exists("scripts/R/utils_enhanced.R")
source_if_exists("Slide-tags/seurat_compat_utils.R")

if (!exists("ensure_dir", mode = "function")) {
  ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_system <- function(cmd) {
  out <- tryCatch(system(cmd, intern = TRUE, ignore.stderr = TRUE), error = function(e) character(0))
  if (length(out) == 0) NA_character_ else paste(out, collapse = "\n")
}

safe_read_rds_summary <- function(path) {
  out <- list(path = path, exists = file.exists(path))
  if (!out$exists) return(out)

  obj <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(obj, "error")) {
    out$read_error <- conditionMessage(obj)
    return(out)
  }

  out$class <- class(obj)[1]
  if (inherits(obj, "Seurat")) {
    out$cells <- ncol(obj)
    out$features <- nrow(obj)
    out$assays <- names(obj@assays)
    out$reductions <- names(obj@reductions)
    out$meta_cols <- colnames(obj@meta.data)
    if ("week" %in% out$meta_cols) {
      out$week_values <- sort(unique(stats::na.omit(obj@meta.data$week)))
    }
    for (cand in c("sample_id", "orig.ident", "donor_id", "biosample_id", "gestational_week")) {
      if (cand %in% out$meta_cols) {
        out[[paste0(cand, "_unique")]] <- head(unique(as.character(obj@meta.data[[cand]])), 20)
      }
    }
  }
  out
}

ensure_dir(DIR_LOGS)
ensure_dir("output_share/logs")

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
report_txt <- file.path("output_share/logs", paste0("pipeline_diagnostic_", stamp, ".txt"))
report_json <- file.path("output_share/logs", paste0("pipeline_diagnostic_", stamp, ".json"))

expected_files <- list(
  reference = file.path(DIR_OBJS, "multiome_reference_processed.rds"),
  slidetags_primary = file.path(DIR_OBJS, "slidetags_mapped.rds"),
  slidetags_compat = file.path(DIR_OBJS, "slidetags_mapped_to_multiome.rds"),
  starmap_primary = file.path(DIR_OBJS, "starmap_mapped.rds"),
  starmap_compat = file.path(DIR_OBJS, "starmap_mapped_to_multiome.rds"),
  starmap_harmonized = file.path(DIR_OBJS, "starmap_harmonized.rds"),
  slidetags_harmonized = file.path(DIR_OBJS, "slidetags_harmonized.rds")
)

script_checks <- c(
  "scripts/02_preprocess/02A_preprocess_multiome_reference.R",
  "scripts/03_mapping/03A_map_slidetags_to_multiome.R",
  "scripts/03_mapping/03B_map_starmap_to_multiome.R",
  "scripts/03_mapping/03C_harmonize_celltype_labels.R",
  "scripts/04_timecourse/04A_gene_of_interest_timecourse.R",
  "scripts/04_timecourse/04A_gene_of_interest_timecourse_ENHANCED.R",
  "scripts/R/utils.R",
  "scripts/R/utils_enhanced.R"
)

diagnostic <- list(
  timestamp = as.character(Sys.time()),
  wd = getwd(),
  r_version = R.version.string,
  seurat_version = if (ok_seurat) as.character(utils::packageVersion("Seurat")) else NA_character_,
  seurat_object_version = if (requireNamespace("SeuratObject", quietly = TRUE)) as.character(utils::packageVersion("SeuratObject")) else NA_character_,
  git_head = safe_system("git rev-parse --short HEAD"),
  git_branch = safe_system("git rev-parse --abbrev-ref HEAD"),
  scripts_present = setNames(file.exists(script_checks), script_checks),
  objects = lapply(expected_files, safe_read_rds_summary),
  key_functions_loaded = c(
    print_seurat_diagnostic = exists("print_seurat_diagnostic", mode = "function"),
    ensure_week_column = exists("ensure_week_column", mode = "function"),
    ensure_spatial_coords = exists("ensure_spatial_coords", mode = "function"),
    log_msg_enhanced = exists("log_msg_enhanced", mode = "function")
  )
)

# Human-readable report
lines <- c(
  "=== PIPELINE DIAGNOSTIC REPORT ===",
  paste("Timestamp:", diagnostic$timestamp),
  paste("Working dir:", diagnostic$wd),
  paste("R:", diagnostic$r_version),
  paste("Seurat:", diagnostic$seurat_version),
  paste("SeuratObject:", diagnostic$seurat_object_version),
  paste("Git branch:", diagnostic$git_branch),
  paste("Git head:", diagnostic$git_head),
  "",
  "-- Key functions loaded --"
)
for (nm in names(diagnostic$key_functions_loaded)) {
  lines <- c(lines, paste(" ", nm, "=", diagnostic$key_functions_loaded[[nm]]))
}

lines <- c(lines, "", "-- Expected script files --")
for (p in names(diagnostic$scripts_present)) {
  lines <- c(lines, paste(" ", p, "=", diagnostic$scripts_present[[p]]))
}

lines <- c(lines, "", "-- Expected object files --")
for (nm in names(diagnostic$objects)) {
  obj <- diagnostic$objects[[nm]]
  lines <- c(lines, paste0("[", nm, "] ", obj$path), paste("  exists:", obj$exists))
  if (!is.null(obj$read_error)) {
    lines <- c(lines, paste("  read_error:", obj$read_error))
  }
  if (!is.null(obj$class)) lines <- c(lines, paste("  class:", obj$class))
  if (!is.null(obj$cells)) lines <- c(lines, paste("  cells:", obj$cells, "features:", obj$features))
  if (!is.null(obj$assays)) lines <- c(lines, paste("  assays:", paste(obj$assays, collapse = ", ")))
  if (!is.null(obj$reductions)) lines <- c(lines, paste("  reductions:", paste(obj$reductions, collapse = ", ")))
  if (!is.null(obj$week_values)) lines <- c(lines, paste("  week values:", paste(obj$week_values, collapse = ", ")))
  for (cand in c("sample_id_unique", "orig.ident_unique", "donor_id_unique", "biosample_id_unique", "gestational_week_unique")) {
    if (!is.null(obj[[cand]])) lines <- c(lines, paste(" ", cand, ":", paste(obj[[cand]], collapse = ", ")))
  }
  lines <- c(lines, "")
}

writeLines(lines, report_txt)
if (ok_jsonlite) {
  jsonlite::write_json(diagnostic, report_json, auto_unbox = TRUE, pretty = TRUE, na = "null")
}

cat("Wrote:\n", report_txt, "\n")
if (ok_jsonlite) cat(report_json, "\n")

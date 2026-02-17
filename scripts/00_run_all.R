#!/usr/bin/env Rscript

# Run pipeline modules in sequence and write reproducibility logs.

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

parse_args <- function(args) {
  out <- list(steps = NULL)
  for (a in args) {
    if (grepl("^--steps=", a)) out$steps <- sub("^--steps=", "", a)
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

step_to_script <- c(
  preprocess = "scripts/02_preprocess/02A_preprocess_multiome_reference.R",
  mapping    = "scripts/03_mapping/03A_map_slidetags_to_multiome.R",
  timecourse = "scripts/04_timecourse/04A_timecourse_analysis.R",
  spatial    = "scripts/05_spatial/05A_spatial_analysis.R",
  publish    = "scripts/99_publish_results.R"
)

default_steps <- names(step_to_script)
steps <- if (is.null(args$steps) || !nzchar(args$steps)) default_steps else strsplit(args$steps, ",", fixed = TRUE)[[1]]
steps <- trimws(steps)

invalid_steps <- setdiff(steps, names(step_to_script))
if (length(invalid_steps) > 0) {
  stop("Unknown step(s): ", paste(invalid_steps, collapse = ", "),
       "\nValid steps: ", paste(names(step_to_script), collapse = ", "))
}

ensure_dir("output/logs")
ensure_dir("output_share/logs")

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path("output/logs", paste0("00_run_all_", run_id, ".log"))
session_file <- file.path("output_share/logs", paste0("sessionInfo_", run_id, ".txt"))

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Pipeline start. Steps: ", paste(steps, collapse = " -> "))

for (s in steps) {
  script <- step_to_script[[s]]
  if (!file.exists(script)) {
    stop("Missing script for step '", s, "': ", script)
  }
  log_msg("Running step '", s, "' via ", script)
  source(script, local = new.env(parent = globalenv()))
  log_msg("Completed step '", s, "'.")
}

writeLines(capture.output(sessionInfo()), session_file)
log_msg("Wrote session info: ", session_file)
log_msg("Pipeline complete.")

#!/usr/bin/env Rscript

# Copy selected outputs into output_share/ and write run manifest.

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

run_cmd <- function(cmd) {
  out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
  if (length(out) == 0) return(NA_character_)
  trimws(out[[1]])
}

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\\\"', x)
  x <- gsub("\n", "\\\\n", x)
  x
}

copy_selected <- function(files, to_dir) {
  if (length(files) == 0) return(character(0))
  ensure_dir(to_dir)
  copied <- character(0)
  for (f in files) {
    target <- file.path(to_dir, basename(f))
    file.copy(f, target, overwrite = TRUE)
    copied <- c(copied, target)
  }
  copied
}

ensure_dir("output_share/figures")
ensure_dir("output_share/tables")
ensure_dir("output_share/logs")

fig_files <- if (dir.exists("output/figures")) {
  list.files("output/figures", pattern = "\\.png$", full.names = TRUE)
} else character(0)

# Keep table publish rules simple and reproducible.
table_patterns <- c("summary", "celltype", "mapping", "timecourse", "spatial")
table_candidates <- if (dir.exists("output/tables")) {
  list.files("output/tables", pattern = "\\.(csv|tsv)$", full.names = TRUE, ignore.case = TRUE)
} else character(0)
if (length(table_candidates) > 0) {
  keep <- Reduce(`|`, lapply(table_patterns, grepl, x = basename(table_candidates), ignore.case = TRUE))
  table_files <- table_candidates[keep]
} else {
  table_files <- character(0)
}

copied_figs <- copy_selected(fig_files, "output_share/figures")
copied_tables <- copy_selected(table_files, "output_share/tables")

manifest_path <- "output_share/logs/run_manifest.json"
manifest <- c(
  "{",
  paste0('  "generated_at": "', format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), '",'),
  paste0('  "git_commit": "', json_escape(run_cmd("git rev-parse --short HEAD")), '",'),
  '  "publish_rules": {',
  '    "figures": "output/figures/*.png",',
  '    "tables": "output/tables/*(summary|celltype|mapping|timecourse|spatial)*.(csv|tsv)"',
  '  },',
  paste0('  "n_figures": ', length(copied_figs), ','),
  paste0('  "n_tables": ', length(copied_tables), ','),
  '  "figures": [',
  if (length(copied_figs) > 0) paste0('    "', json_escape(copied_figs), '"', collapse = ",\n") else "",
  '  ],',
  '  "tables": [',
  if (length(copied_tables) > 0) paste0('    "', json_escape(copied_tables), '"', collapse = ",\n") else "",
  '  ]',
  "}"
)
writeLines(manifest, manifest_path)

cat("Published ", length(copied_figs), " figure(s) and ", length(copied_tables), " table(s).\n", sep = "")
cat("Manifest: ", manifest_path, "\n", sep = "")

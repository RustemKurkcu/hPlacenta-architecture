# Core utility fallbacks for pipeline scripts.
# This file is intentionally lightweight so updated script bundles can run
# even if a larger historical utils.R is missing pieces.

`%||%` <- function(a, b) if (!is.null(a)) a else b

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

log_msg <- function(msg, logfile = NULL) {
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  if (!is.null(logfile)) cat(line, "\n", file = logfile, append = TRUE)
  invisible(line)
}

pick_first_present <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

ensure_week_column <- function(obj, candidates = c("week", "gestational_week", "donor_id", "biosample_id", "sample_id")) {
  if (!inherits(obj, "Seurat")) return(obj)
  md <- obj@meta.data
  if ("week" %in% colnames(md) && any(!is.na(md$week))) return(obj)

  src <- pick_first_present(md, candidates)
  if (is.null(src)) {
    obj$week <- NA_real_
    return(obj)
  }

  vals <- as.character(md[[src]])
  wk <- suppressWarnings(as.numeric(vals))
  if (all(is.na(wk))) {
    wk <- suppressWarnings(as.numeric(sub(".*?(\\d+).*", "\\1", vals)))
  }
  obj$week <- wk
  obj
}

ensure_spatial_coords <- function(obj, x_candidates, y_candidates) {
  if (!inherits(obj, "Seurat")) return(obj)
  md <- obj@meta.data
  xcol <- pick_first_present(md, x_candidates)
  ycol <- pick_first_present(md, y_candidates)
  obj$spatial_x_use <- if (!is.null(xcol)) suppressWarnings(as.numeric(md[[xcol]])) else NA_real_
  obj$spatial_y_use <- if (!is.null(ycol)) suppressWarnings(as.numeric(md[[ycol]])) else NA_real_
  obj
}

print_seurat_diagnostic <- function(obj, label = "Seurat object") {
  cat(strrep("=", 70), "\n")
  cat(label, "Diagnostic\n")
  cat(strrep("=", 70), "\n")
  cat("Cells:", ncol(obj), "\n")
  cat("Assays:", paste(names(obj@assays), collapse = ", "), "\n")
  cat("Reductions:", paste(names(obj@reductions), collapse = ", "), "\n")
  cat("Metadata columns:", ncol(obj@meta.data), "\n")
  if ("week" %in% colnames(obj@meta.data)) {
    w <- sort(unique(stats::na.omit(obj@meta.data$week)))
    if (length(w) > 0) cat("Weeks:", paste(w, collapse = ", "), "\n")
  }
  cat(strrep("=", 70), "\n\n")
  invisible(NULL)
}

save_plot <- function(plot_obj, path, w = 8, h = 6, dpi = 300) {
  ggplot2::ggsave(filename = path, plot = plot_obj, width = w, height = h, dpi = dpi)
}

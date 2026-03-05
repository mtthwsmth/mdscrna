#' Order numeric cluster labels in Seurat metadata
#'
#' Reorders factor levels for metadata columns whose names start with any of the supplied
#' \code{prefix} values so that numeric cluster labels (e.g. \code{"0"}, \code{"1"}, \code{"2"})
#' are in increasing order. Errors if any matching column contains non-numeric labels
#' (e.g. \code{"2_1"}).
#'
#' @param seurat_obj A Seurat object.
#' @param prefix Character vector of prefix strings used to match metadata columns.
#'   Default \code{c("SCT_snn_res.", "RNA_snn_res.")}.
#' @param update_ident Optional metadata column name to set as the active identity (Idents).
#'
#' @return The Seurat object with reordered factor levels (and optionally updated Idents).
#' @export
mds_so_order_clusters_numeric <- function(seurat_obj,
                                          prefix = c("SCT_snn_res.", "RNA_snn_res."),
                                          update_ident = NULL) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.", call. = FALSE)
  }

  if (!is.character(prefix) || length(prefix) < 1) {
    stop("prefix must be a non-empty character vector.", call. = FALSE)
  }

  md_cols <- colnames(seurat_obj@meta.data)

  # Match any prefix at start of column name
  pref_rx <- paste0("^(", paste(vapply(prefix, function(p) gsub("([.^$|()*+?{\\[\\]\\\\])", "\\\\\\1", p),
                                       character(1)), collapse = "|"), ")")
  cluster_cols <- grep(pref_rx, md_cols, value = TRUE)

  if (length(cluster_cols) == 0) {
    warning(
      paste0("No columns found matching prefixes: ", paste(prefix, collapse = ", ")),
      call. = FALSE
    )
    return(seurat_obj)
  }

  for (col in cluster_cols) {
    current_values <- seurat_obj@meta.data[[col]]

    if (is.factor(current_values) || is.character(current_values)) {
      vals <- as.character(current_values)

      suppressWarnings(num <- as.numeric(vals))

      if (any(is.na(num) & !is.na(vals))) {
        bad <- unique(vals[is.na(num) & !is.na(vals)])
        stop(sprintf(
          "Column '%s' contains non-numeric cluster labels (examples: %s).",
          col, paste(utils::head(bad, 5), collapse = ", ")
        ), call. = FALSE)
      }

      level_chars <- as.character(sort(unique(num)))
      seurat_obj@meta.data[[col]] <- factor(vals, levels = level_chars)

      message(paste0("Reordered column: ", col))
    } else {
      message(paste0("Skipping column '", col, "' - not a factor or character"))
    }
  }

  if (!is.null(update_ident)) {
    if (update_ident %in% colnames(seurat_obj@meta.data)) {
      Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[update_ident]]
      message(paste0("Updated active identity to: ", update_ident))
    } else {
      warning(
        paste0("Column '", update_ident, "' not found in metadata. Active identity not changed."),
        call. = FALSE
      )
    }
  }

  seurat_obj
}

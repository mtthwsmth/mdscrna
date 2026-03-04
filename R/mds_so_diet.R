#' Clean a Seurat object down to RNA counts + (filtered) meta.data
#'
#' Creates a new Seurat object containing only the specified assay counts layer and a cleaned
#' meta.data table (optionally removing clustering columns).
#'
#' If the requested layer (default \code{"counts"}) is not present but split layers are detected
#' (e.g. \code{counts.sample1}, \code{counts.sample2}), the function runs
#' \code{SeuratObject::JoinLayers()} internally to create the requested layer and prints a message.
#'
#' @param so A Seurat object.
#' @param assay Assay to extract counts from. Default \code{"RNA"}.
#' @param layer Which layer to treat as counts. Default \code{"counts"}.
#' @param remove_clustering Logical; remove clustering columns from meta.data. Default \code{TRUE}.
#' @param clustering_prefixes Character vector of clustering prefixes to remove. Default
#'   \code{c("SCT_snn_res", "RNA_snn_res")}.
#' @param additional_remove_patterns Optional character vector of additional regex patterns.
#' @param verbose Logical; emit messages. Default \code{TRUE}.
#'
#' @return A new Seurat object with counts and cleaned meta.data.
#' @export
mds_so_diet <- function(
    so,
    assay = "RNA",
    layer = "counts",
    remove_clustering = TRUE,
    clustering_prefixes = c("SCT_snn_res", "RNA_snn_res"),
    additional_remove_patterns = NULL,
    verbose = TRUE
) {
  if (!inherits(so, "Seurat")) stop("'so' must be a Seurat object.")
  if (!assay %in% names(so@assays)) stop(sprintf("Assay '%s' not found.", assay))

  md <- so@meta.data

  remove_patterns <- character(0)
  if (isTRUE(remove_clustering) && length(clustering_prefixes) > 0) {
    remove_patterns <- c(remove_patterns, paste0("^(", paste(clustering_prefixes, collapse = "|"), ")\\."))
  }
  if (!is.null(additional_remove_patterns) && length(additional_remove_patterns) > 0) {
    remove_patterns <- c(remove_patterns, additional_remove_patterns)
  }

  if (length(remove_patterns) > 0 && ncol(md) > 0) {
    to_remove <- Reduce(
      `|`,
      lapply(remove_patterns, function(p) grepl(p, colnames(md))),
      init = rep(FALSE, ncol(md))
    )
    n_removed <- sum(to_remove)
    if (n_removed > 0) {
      if (isTRUE(verbose)) {
        message(sprintf("Removing %d meta.data column(s) matching clustering/pattern rules.", n_removed))
      }
      md <- md[, !to_remove, drop = FALSE]
    }
  }

  layers <- SeuratObject::Layers(so[[assay]])

  if (!(layer %in% layers)) {
    split_layers <- grep(paste0("^", layer, "\\."), layers, value = TRUE)
    if (length(split_layers) > 0) {
      if (isTRUE(verbose)) {
        message(sprintf(
          "Layer '%s' not found in assay '%s', but split layers detected (%s). Running SeuratObject::JoinLayers()...",
          layer, assay, paste(head(split_layers, 5), collapse = ", ")
        ))
      }
      so <- SeuratObject::JoinLayers(so, assay = assay)
      layers <- SeuratObject::Layers(so[[assay]])
    }
  }

  if (!(layer %in% layers)) {
    stop(sprintf(
      "Layer '%s' not found in assay '%s'. Available layers: %s",
      layer, assay, paste(layers, collapse = ", ")
    ))
  }

  counts <- SeuratObject::LayerData(so, assay = assay, layer = layer)
  new <- Seurat::CreateSeuratObject(counts = counts, assay = assay)

  # ---- robust alignment check ----
  idx <- match(colnames(new), rownames(md))
  if (anyNA(idx)) {
    missing <- colnames(new)[is.na(idx)]
    stop(sprintf(
      "Cell name mismatch between counts matrix and meta.data: %d cell(s) missing (examples: %s).",
      length(missing),
      paste(utils::head(missing, 5), collapse = ", ")
    ))
  }
  md <- md[idx, , drop = FALSE]

  new <- Seurat::AddMetaData(new, metadata = md)
  new
}

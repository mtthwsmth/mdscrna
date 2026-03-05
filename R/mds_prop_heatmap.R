#' Proportion heatmap from metadata
#'
#' Plots the row-wise proportions (percent) of one or two metadata fields within groups
#' (e.g., clusters) for a Seurat or SingleCellExperiment object.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param features Character vector of length 1 or 2 giving metadata column name(s).
#' @param group.by Metadata column name to group rows by. Default "seurat_clusters".
#' @param scale Logical; if TRUE, z-score values before plotting. Default FALSE.
#' @param colors Two colors defining the fill gradient. Default c("white", "red").
#' @param title Plot title. Default "Proportion Heatmap".
#'
#' @return A ggplot object. The underlying matrix is available via
#'   \code{attr(plot, "prop_table")}.
#'
#' @export
mds_prop_heatmap <- function(
    object,
    features,
    group.by = "seurat_clusters",
    scale = FALSE,
    colors = c("white", "red"),
    title = "Proportion Heatmap"
) {
  # deps
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for mds_prop_heatmap.")
  }

  # metadata: support Seurat OR SingleCellExperiment
  meta_data <- NULL
  if (inherits(object, "Seurat")) {
    meta_data <- object@meta.data
  } else if (inherits(object, "SingleCellExperiment")) {
    meta_data <- as.data.frame(SummarizedExperiment::colData(object))
  } else {
    stop("'object' must be a Seurat or SingleCellExperiment object.")
  }

  # checks
  missing_features <- features[!features %in% colnames(meta_data)]
  if (length(missing_features) > 0) {
    stop("The following features are not found in metadata: ",
         paste(missing_features, collapse = ", "))
  }
  if (!group.by %in% colnames(meta_data)) {
    stop("Grouping variable '", group.by, "' not found in metadata")
  }

  # build table (include NA)
  if (length(features) == 1) {
    col_vec <- meta_data[[features[1]]]
  } else {
    # match previous behavior: combine first two features with "_"
    col_vec <- paste(meta_data[[features[1]]],
                     meta_data[[features[2]]], sep = "_")
  }
  grp_vec <- meta_data[[group.by]]

  prop_table <- table(grp_vec, col_vec, useNA = "ifany")

  # proportions by row (each group sums to 100%)
  rs <- rowSums(prop_table)
  prop_mat <- suppressWarnings(sweep(prop_table, 1, rs, "/") * 100)
  prop_mat[is.nan(prop_mat)] <- NA_real_

  # optional column-wise scaling (mimics prior `scale(prop_table)`)
  if (isTRUE(scale)) {
    z <- scale(prop_mat)
    dimnames(z) <- dimnames(prop_mat)
    prop_mat <- z
  }

  # tidy for ggplot
  df <- as.data.frame(as.table(prop_mat), stringsAsFactors = FALSE)
  colnames(df) <- c("group", "label", "value")

  # keep ordering stable & readable
  df$group <- factor(df$group, levels = rownames(prop_mat))
  df$label <- factor(df$label, levels = colnames(prop_mat))

  # text labels (match previous %.1f)
  df$lab <- ifelse(is.na(df$value), "", sprintf("%.1f", df$value))

  # colors
  pal <- grDevices::colorRampPalette(colors)(100)

  # build plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = label, y = group, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = lab), size = 3) +
    ggplot2::scale_fill_gradientn(colors = pal, na.value = "grey90") +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = if (scale) "scaled" else "%") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  # attach underlying matrix
  attr(p, "prop_table") <- prop_mat
  p
}

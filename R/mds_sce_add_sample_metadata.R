#' Add sample-level metadata to a SingleCellExperiment
#'
#' Left-joins a sample-level \code{sample_meta} table onto \code{colData(sce)}. This is intended
#' for the common case where \code{sample_meta} has one row per sample (unique sample IDs),
#' and \code{colData(sce)} has many rows per sample (cells), so the sample metadata is
#' replicated across cells from the same sample.
#'
#' The function hard-checks that:
#' \itemize{
#'   \item \code{sample_meta[[metadata_by]]} contains unique keys (to prevent row duplication).
#'   \item Neither \code{colData(sce)[[coldata_by]]} nor \code{sample_meta[[metadata_by]]}
#'         may contain \code{NA} values in the join key.
#'   \item No columns in \code{sample_meta} (other than the join key) share a name with
#'         existing \code{colData(sce)} columns — this would introduce ambiguous
#'         \code{.x}/\code{.y} suffixed columns and is treated as an error.
#'   \item When \code{coldata_by != metadata_by}, neither key column name may appear in
#'         the opposite data frame (which would also produce \code{.x}/\code{.y} suffixes).
#'   \item Cell order is preserved (validated via a temporary cell-id column).
#' }
#'
#' @param sce A \code{SingleCellExperiment}.
#' @param sample_meta A \code{data.frame} (or coercible to one) containing sample-level
#'   metadata to add.
#' @param coldata_by Column name in \code{colData(sce)} to join on. Default is \code{"Sample"}.
#' @param metadata_by Column name in \code{sample_meta} to join on. Default is \code{"Sample"}.
#' @param keep_key_cols Logical; if \code{FALSE} and \code{coldata_by != metadata_by},
#'   drops the join-key column coming from \code{sample_meta} after the join (the
#'   \code{colData} key is retained). Default \code{TRUE}.
#'
#' @return The input \code{sce} with updated \code{colData}.
#'
#' @importFrom methods is
#' @export
mds_sce_add_sample_metadata <- function(
    sce,
    sample_meta,
    coldata_by = "Sample",
    metadata_by = "Sample",
    keep_key_cols = TRUE
) {
  if (!methods::is(sce, "SingleCellExperiment")) {
    stop("'sce' must be a SingleCellExperiment object.")
  }
  if (!is.data.frame(sample_meta)) sample_meta <- as.data.frame(sample_meta)

  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  cell_names <- colnames(sce)

  # Ensure existing colData is aligned to the SCE columns (cells)
  rn <- rownames(cd)
  if (is.null(rn) || anyNA(rn) || !all(rn == cell_names)) {
    stop("colData(sce) rownames must be present and identical to colnames(sce) (in the same order).")
  }

  if (!coldata_by %in% names(cd)) {
    stop(sprintf("Join key '%s' not found in colData(sce).", coldata_by))
  }
  if (!metadata_by %in% names(sample_meta)) {
    stop(sprintf("Join key '%s' not found in sample_meta.", metadata_by))
  }

  # Hard check: NA keys are not permitted on either side
  if (anyNA(cd[[coldata_by]])) {
    stop(sprintf(
      "Join key '%s' in colData(sce) contains NA values. Remove or resolve these before joining.",
      coldata_by
    ))
  }
  if (anyNA(sample_meta[[metadata_by]])) {
    stop(sprintf(
      "Join key '%s' in sample_meta contains NA values. Remove or resolve these before joining.",
      metadata_by
    ))
  }

  # Hard check: sample_meta keys must be unique (sample-level table)
  key <- sample_meta[[metadata_by]]
  if (anyDuplicated(key)) {
    dup <- unique(key[duplicated(key)])
    stop(sprintf(
      "sample_meta must have unique '%s' values; found %d duplicated key value(s) (examples: %s).",
      metadata_by,
      length(dup),
      paste(utils::head(dup, 5), collapse = ", ")
    ))
  }

  # Hard check: no incoming columns may share a name with existing colData columns
  existing_cols <- setdiff(names(cd), coldata_by)
  incoming_cols <- setdiff(names(sample_meta), metadata_by)
  conflicts <- intersect(existing_cols, incoming_cols)
  if (length(conflicts) > 0) {
    stop(sprintf(
      "The following columns in sample_meta already exist in colData(sce): %s. Remove or rename them before joining.",
      paste(conflicts, collapse = ", ")
    ))
  }

  # Hard check: when key names differ, neither key column name may appear in the opposite data frame
  if (coldata_by != metadata_by) {
    if (metadata_by %in% names(cd)) {
      stop(sprintf(
        "metadata_by ('%s') already exists as a column in colData(sce). Rename it before joining to avoid .x/.y suffixes.",
        metadata_by
      ))
    }
    if (coldata_by %in% names(sample_meta)) {
      stop(sprintf(
        "coldata_by ('%s') already exists as a column in sample_meta. Rename it before joining to avoid .x/.y suffixes.",
        coldata_by
      ))
    }
  }

  # Preserve cell identity across dplyr join (since rownames are not preserved)
  .mds_cell_id_col <- ".mds_cell_id"
  if (.mds_cell_id_col %in% names(cd) || .mds_cell_id_col %in% names(sample_meta)) {
    stop(sprintf("Temporary column '%s' already exists; please rename/remove it and retry.", .mds_cell_id_col))
  }
  cd[[.mds_cell_id_col]] <- cell_names

  by <- stats::setNames(metadata_by, coldata_by)
  cd2 <- dplyr::left_join(cd, sample_meta, by = by)

  if (nrow(cd2) != nrow(cd)) {
    stop(sprintf(
      "Join unexpectedly changed row count (expected %d, got %d).",
      nrow(cd), nrow(cd2)
    ))
  }

  # Validate order via the temporary cell-id column
  if (!.mds_cell_id_col %in% names(cd2)) {
    stop("Internal error: temporary cell-id column missing after join.")
  }
  if (!all(cd2[[.mds_cell_id_col]] == cell_names)) {
    stop("Cell order mismatch after join: joined colData no longer matches colnames(sce).")
  }

  # Restore rownames deterministically and drop temporary column
  rownames(cd2) <- cd2[[.mds_cell_id_col]]
  cd2[[.mds_cell_id_col]] <- NULL

  if (!keep_key_cols && coldata_by != metadata_by) {
    cd2[[metadata_by]] <- NULL
  }

  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cd2)
  sce
}

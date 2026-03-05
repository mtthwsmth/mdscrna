#' Union variable features across groups for Harmony-style workflows
#'
#' Computes variable features separately within each group defined by \code{split_by},
#' takes the union across groups, and either assigns them to the Seurat object
#' (\code{VariableFeatures(so) <- union_features}) or returns the feature vector.
#'
#' This follows the approach used in the Harmony Seurat v5 guide: split cells by a
#' metadata field, run \code{FindVariableFeatures()} per group, take \code{unique(unlist(...))}
#' and assign the union back to \code{VariableFeatures(object)} before downstream steps
#' like scaling, PCA, and \code{RunHarmony()}.  [oai_citation:0‡immunogenomics.r-universe.dev](https://immunogenomics.r-universe.dev/harmony/doc/Seurat.html)
#'
#'#' @details
#' This helper does not run normalization. For typical workflows, run
#' \code{Seurat::NormalizeData()} (and any required preprocessing) before calling this
#' function so \code{FindVariableFeatures()} operates on the intended data.
#'
#' @param so A Seurat object.
#' @param split_by Metadata column name used to split cells into groups.
#' @param nfeatures Number of variable features to select per group. Default 3000.
#' @param selection.method Variable feature selection method passed to
#'   \code{Seurat::FindVariableFeatures()}. Default \code{"vst"}.
#' @param assay Optional assay name to set as the default assay before selecting features.
#'   If \code{NULL}, uses the current default assay.
#' @param ret Return mode: \code{"object"} returns the Seurat object with union features set;
#'   \code{"features"} returns the union feature vector. Default \code{"object"}.
#' @param verbose Logical; print messages. Default \code{TRUE}.
#'
#' @return If \code{ret = "object"}, a Seurat object with \code{VariableFeatures(so)}
#'   set to the union of group-wise variable features. If \code{ret = "features"},
#'   a character vector of feature names.
#'
#' @seealso \code{\link[Seurat:FindVariableFeatures]{Seurat::FindVariableFeatures}},
#'   \code{\link[harmony:RunHarmony]{harmony::RunHarmony}}
#'
#' @export
mds_so_set_union_vf <- function(
    so,
    split_by,
    nfeatures = 3000,
    selection.method = "vst",
    assay = NULL,
    ret = "object",
    verbose = TRUE
) {
  if (!inherits(so, "Seurat")) {
    stop("Argument 'so' must be a Seurat object.")
  }
  if (!split_by %in% colnames(so@meta.data)) {
    stop(sprintf("Column '%s' not in meta.data.", split_by))
  }
  if (!ret %in% c("object", "features")) {
    stop("ret must be 'object' or 'features'.")
  }
  if (!is.null(assay)) Seurat::DefaultAssay(so) <- assay

  keep <- !is.na(so@meta.data[[split_by]])
  if (!all(keep)) {
    if (isTRUE(verbose)) {
      message(sprintf(
        "Dropping %d cell(s) with NA in '%s'.",
        sum(!keep), split_by
      ))
    }
    so <- so[, keep]
  }

  groups <- split(rownames(so@meta.data), so@meta.data[[split_by]])
  if (isTRUE(verbose)) {
    message(sprintf("Found %d group(s).", length(groups)))
  }

  vf_list <- lapply(groups, function(cells_use) {
    x <- so[, cells_use]
    if (ncol(x) < 50 && isTRUE(verbose)) {
      message(sprintf(
        "Group has %d cells; VST features may be noisy.",
        ncol(x)
      ))
    }
    x <- Seurat::FindVariableFeatures(
      x,
      selection.method = selection.method,
      nfeatures = nfeatures,
      verbose = FALSE
    )
    Seurat::VariableFeatures(x)
  })

  vfs <- unique(unlist(vf_list, use.names = FALSE))
  if (length(vfs) == 0) {
    stop("No variable features were identified.")
  }

  if (isTRUE(verbose)) {
    message(sprintf("Collected %d union feature(s).", length(vfs)))
  }

  if (ret == "features") return(vfs)

  Seurat::VariableFeatures(so) <- vfs
  if (isTRUE(verbose)) {
    message("Assigned union variable features to the object.")
  }
  so
}

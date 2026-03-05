#' Run GO enrichment on up/down gene sets
#'
#' Runs \code{clusterProfiler::enrichGO()} on the \code{up} and \code{down} results in a
#' DGE split list (e.g. output of \code{mds_dge_split}). If any terms are enriched,
#' the result is augmented with \code{enrichplot::pairwise_termsim()}.
#'
#' @param df.ls A named list containing \code{up} and \code{down} data.frames. Row names
#'   should be gene identifiers (SYMBOL by default).
#' @param ontology GO ontology to use. One of \code{"BP"}, \code{"MF"}, or \code{"CC"}.
#'   Default \code{"BP"}.
#' @param universe Character vector of gene IDs defining the background universe.
#' @param OrgDb An \code{OrgDb} object (e.g. \code{org.Hs.eg.db::org.Hs.eg.db}).
#' @param ... Additional arguments passed to \code{clusterProfiler::enrichGO()}.
#'
#' @return A named list containing \code{ego.up} and/or \code{ego.down}. If the corresponding
#'   input set is empty, that element is omitted.
#'
#' @export
mds_enrichgo_up_down <- function(df.ls, ontology = "BP", universe, OrgDb, ...) {
  ego.up <- NULL
  ego.down <- NULL

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required.")
  }
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required for pairwise_termsim().")
  }

  if (nrow(df.ls$up) > 0) {
    ego.up <- clusterProfiler::enrichGO(
      rownames(df.ls$up),
      keyType = "SYMBOL",
      OrgDb = OrgDb,
      ont = ontology,
      universe = universe,
      ...
    )
    if (nrow(as.data.frame(ego.up)) > 0) {
      ego.up <- enrichplot::pairwise_termsim(ego.up)
    }
  }

  if (nrow(df.ls$down) > 0) {
    ego.down <- clusterProfiler::enrichGO(
      rownames(df.ls$down),
      keyType = "SYMBOL",
      OrgDb = OrgDb,
      ont = ontology,
      universe = universe,
      ...
    )
    if (nrow(as.data.frame(ego.down)) > 0) {
      ego.down <- enrichplot::pairwise_termsim(ego.down)
    }
  }

  result <- list()
  if (!is.null(ego.up)) result$ego.up <- ego.up
  if (!is.null(ego.down)) result$ego.down <- ego.down
  result
}

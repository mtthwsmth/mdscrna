#' Split DGE results into up- and down-regulated sets
#'
#' Given a DESeq2 or edgeR results table, filters by fold-change and adjusted
#' p-value and returns the full filtered table plus up/down subsets.
#' Fold-change cutoff is specified in linear space and converted to log2.
#' Rows with NA values are removed before filtering.
#'
#' @param result A data.frame-like object with one row per gene.
#' @param foldChangeName Column name containing log2 fold-change.
#' @param foldChangeCutoff Fold-change cutoff in linear space (e.g. 1.5).
#' @param pvalName Column name containing adjusted p-values (e.g. "padj" or "FDR").
#' @param pvalCutoff Adjusted p-value cutoff (e.g. 0.05).
#' @param expressionCutoffName Column name containing an expression metric (e.g. "logCPM").
#' @param expressionCutoff Optional numeric cutoff for the expression metric; if NULL, no expression filtering is applied.
#'
#' @return A list with three elements: \code{res}, \code{up}, and \code{down}.
#' @export

mds_dge_split <- function(result, foldChangeName = "logFC", foldChangeCutoff = (1.5), pvalName = "FDR", pvalCutoff = 0.05, expressionCutoffName = "logCPM", expressionCutoff = NULL) {
  # subset the results to only include genes that pass the fold change and p value cutoffs
  x <- na.omit(result) #both DESeq2 and edgeR downstream of pseudoBulkDGE return resutls with NAs in them, so get rid of those.
  cutoff <- log2(foldChangeCutoff)
  # if an expression cutoff is set, use it, otherwise don't filter on expression and just on logFC and pvalue
  if(!is.null(expressionCutoff)) {
    up <- x[x[[foldChangeName]] > cutoff & x[[pvalName]] < pvalCutoff & x[[expressionCutoffName]] > expressionCutoff,]
    down <- x[x[[foldChangeName]] < -cutoff & x[[pvalName]] < pvalCutoff & x[[expressionCutoffName]] > expressionCutoff,]
  } else {
    up <- x[x[[foldChangeName]] > cutoff & x[[pvalName]] < pvalCutoff,]
    down <- x[x[[foldChangeName]] < -cutoff & x[[pvalName]] < pvalCutoff,]
  }
  return(list(res = x, up = up, down = down))
}

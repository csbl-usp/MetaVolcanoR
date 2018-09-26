#' A DEG definition function
#'
#' This function creates a new variable indicating DEGs as -1, 0, 1 based on the
#' user-defined foldchange and p-value criteria
#' @param geo2r data.frame/data.table with DE results where lines are genes
#' @param columnStatistics column name of the pvalue variable <strings>
#' @param columnlogFC column name of the foldchange variable <string>
#' @param pvalue thresholds <ex: 0.05>
#' @param logfc foldchange threshold <ex: 0.4>
#' @keywords DEG definition
#' @export
#' @examples
#' deg.def()
deg.def <- function(geo2r, columnStatistics, columnlogFC, pvalue, logfc) {

	mutate(geo2r, deg = ifelse(as.numeric(!!as.name(columnStatistics)) < pvalue & as.numeric(!!as.name(columnlogFC)) < (-1*logfc), -1,
				   ifelse(as.numeric(!!as.name(columnStatistics)) < pvalue & as.numeric(!!as.name(columnlogFC)) > logfc, 1, 0)))

}

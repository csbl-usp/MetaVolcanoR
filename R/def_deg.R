#' A DEG definition function
#'
#' This function creates a new variable indicating DEGs as -1, 0, 1 based on the
#' user-defined foldchange and p-value criteria
#' @param diffexp data.frame/data.table with DE results where lines are genes
#' @param pcriteria column name of the pvalue variable <strings>
#' @param foldchangecol column name of the foldchange variable <string>
#' @param pv thresholds <ex: 0.05>
#' @param fc foldchange threshold <ex: 0.4>
#' @keywords DEG definition
#' @export
#' @examples
#' deg.def()
deg.def <- function(diffexp, pcriteria, foldchangecol, pv, fc) {
	if(pcriteria != "pv" && foldchangecol != "fc") {
		dplyr::mutate(diffexp, deg = ifelse(as.numeric(!!as.name(pcriteria)) < pv & as.numeric(!!as.name(foldchangecol)) < (-1*fc), -1,
				   ifelse(as.numeric(!!as.name(pcriteria)) < pv & as.numeric(!!as.name(foldchangecol)) > fc, 1, 0)))
	} else {
		stop("Oops! Sorry, could you please rename your pcriteria and/or foldchangecol variables?")
	}
}

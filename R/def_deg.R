#' A DEG definition function
#'
#' This function creates a new variable indicating DEGs as -1, 0, 1 based on the
#' user-defined foldchange and p-value criteria
#' @param diffexp data.frame/data.table with differential expression results
#' @param pcriteria column name of the pvalue variable <strings>
#' @param foldchangecol column name of the foldchange variable <string>
#' @param pv pvalue threshold <double>
#' @param fc foldchange threshold <double>
#' @keywords DEG definition
#' @export
#' @examples
#' deg_def()
deg_def <- function(diffexp, pcriteria, foldchangecol, pv, fc) {
    if(pcriteria != "pv" && foldchangecol != "fc") {
        dplyr::mutate(diffexp, 
	    deg = ifelse(as.numeric(!!as.name(pcriteria)) < pv 
			 & as.numeric(!!as.name(foldchangecol)) < (-1*fc), -1,
	          ifelse(as.numeric(!!as.name(pcriteria)) < pv 
			 & as.numeric(!!as.name(foldchangecol)) > fc, 1, 0)))
	} else {
		stop("Oops! Sorry, could you please rename your pcriteria 
		     and/or foldchangecol variables?")
	}
}

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
#' @return \code{data.table/data.frame} with a new deg variable
#' @export
#' @examples
#' data(diffexplist)
#' diffexp <- deg_def(diffexplist[[1]], "pvalue", "Log2FC", 0.05, 0)
#' table(diffexp[['deg']])
deg_def <- function(diffexp, pcriteria, foldchangecol, pv, fc) {

    if(pcriteria == "pv") {
	
        stop("Oops! Sorry, could you please rename your pcriteria
	     variable?")

    } else if(foldchangecol == "fc") {
	
	 stop("Oops! Sorry, could you please rename your foldchangecol
	      variables")

    } else {

        dplyr::mutate(diffexp, 
	    deg = ifelse(as.numeric(!!as.name(pcriteria)) < pv 
			 & as.numeric(!!as.name(foldchangecol)) < (-1*fc), -1,
	          ifelse(as.numeric(!!as.name(pcriteria)) < pv 
			 & as.numeric(!!as.name(foldchangecol)) > fc, 1, 0)))
    }
}

#' A function to calculate variance from confidence interval limits
#'
#' This function takes the limits of a confidence interval (95%)
#' a calculate a variance
#' @param diffexp data.frame/data.table containing differential expression
#'        results
#' @param llcol column name of the fold change coinfidence interval left limit
#'        name <string>
#' @param rlcol column name of the fold change coinfidence interval left limit
#'        name <string>
#' @param vcol column name of a pre-computed variance column <string> (optional)
#' @param foldchangecol column name of the fold change <string> (used as fallback)
#' @param statcol column name of the test statistic <string> (used as fallback)
#' @keywords variance from confidence interval
#' @return \code{data.table/data.frame} with a new \code{vi} variable
#' @export
#' @examples
#' data(diffexplist)
#' diffexp <- calc_vi(diffexplist[[1]], "CI.L", "CI.R")
#' head(diffexp, 3)
calc_vi <- function(diffexp, llcol=NULL, rlcol=NULL, vcol=NULL,
                    foldchangecol="log2FC", statcol="stat") {
    if (!is.null(llcol) && !is.null(rlcol)) {
        diffexp[["vi"]] <- ((as.numeric(diffexp[[rlcol]]) - as.numeric(diffexp[[llcol]]))/3.92)^2
    } else if (!is.null(vcol)) {
        diffexp[["vi"]] <- as.numeric(diffexp[[vcol]])
    } else {
        diffexp[["vi"]] <- (as.numeric(diffexp[[foldchangecol]]) / as.numeric(diffexp[[statcol]]))^2
    }
    diffexp
}

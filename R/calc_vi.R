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
#' @keywords variance from confidence interval
#' @export
#' @examples
#' calc_vi()
calc_vi <- function(diffexp, llcol, rlcol) {
    diffexp[['vi']] <- apply(diffexp, 1, function(gene) {
        ((as.numeric(gene[rlcol]) - as.numeric(gene[llcol]))/3.92)^2
    })
  diffexp
}

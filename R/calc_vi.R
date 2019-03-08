#' A function to calculate variance from confidence interval limits 
#'
#' This function takes the limits of a confidence interval a calculate a variance
#' @param geo2r data.frame/data.table containing all the inputed GEO2R outputs
#' @keywords variance from confidence interval
#' @export
#' @examples
#' calc.vi()
calc.vi <- function(geo2r) {
  print(dim(geo2r))
  geo2r[['vi']] <- apply(geo2r, 1, function(gene) {
    ((as.numeric(gene['CI.R']) - as.numeric(gene['CI.L']))/3.92)^2 
  })
  geo2r
}
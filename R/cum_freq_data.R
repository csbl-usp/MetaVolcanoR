#' A data formating function for inverse-cummulative DEG distribution
#'
#' This function counts how many genes consistly appears as DE along the inputed GEO2R outputs
#' @param meta_geo2r data.frame/data.table containing all the inputed GEO2R outputs  
#' @param nstud the number of inputed GEO2R outputs  <integer>
#' @keywords format data inverse-cummulative DEG distribution
#' @export
#' @examples
#' cum.freq.data()
cum.freq.data <- function(meta_geo2r, nstud) {
  data.frame(DEGs = sapply(0:nstud, function(idx) {
    length(which(meta_geo2r[['ndeg']] >= idx))
  }),
  ndatasets = 0:nstud
  )
}
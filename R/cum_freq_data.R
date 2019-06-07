#' A data formating function for inverse-cummulative DEG distribution
#'
#' This function counts how many genes consistly appears as DE along the 
#' input studies
#' @param meta_diffexp data.frame/data.table containing all the input studies
#' @param nstud the number of inputed GEO2R outputs  <integer>
#' @keywords format data inverse-cummulative DEG distribution
#' @return \code{data.frame} inverse cummulative distribution
#' @export
#' @examples
#' cum_freq_data()
cum_freq_data <- function(meta_diffexp, nstud) {
    data.frame(DEGs = sapply(0:nstud, function(idx) {
                        length(which(meta_diffexp[['ndeg']] >= idx))
                      }),
               ndatasets = 0:nstud
    )
}

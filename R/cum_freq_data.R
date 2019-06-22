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
#' library(dplyr)
#' data(diffexplist)
#' diffexp <- lapply(diffexplist, function(...) deg_def(..., "pvalue", 
#'            "Log2FC", 0.05, 0))
#' diffexp <- rename_col(diffexp, "Symbol")
#' meta_diffexp <- Reduce(function(...) merge(..., by = "Symbol", all = TRUE),
#'            diffexp)
#' meta_diffexp %>%
#' dplyr::select(dplyr::matches("deg_")) %>%
#'     data.matrix -> n_deg
#' meta_diffexp[['ndeg']] <- rowSums(n_deg^2, na.rm = TRUE)
#' cfd <- cum_freq_data(meta_diffexp, length(diffexplist))
#' head(cfd, 3)
cum_freq_data <- function(meta_diffexp, nstud) {
    data.frame("DEGs" = vapply(0:nstud, function(idx) {
                        length(which(meta_diffexp[['ndeg']] >= idx))
                      }, numeric(1)),
               "ndatasets" = 0:nstud
    )
}

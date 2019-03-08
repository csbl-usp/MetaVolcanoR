#' A function setting data format for DEG barplot visualization
#'
#' This function summarize the variable deg from the deg.def() function to visualize as barplots the number of DEGs across the inputed GEO2R outputs
#' @param geo2r_res list of data.frame/data.table (s) output of the deg.def() function <list>
#' @keywords set DEG barplot data format
#' @export
#' @examples
#' set.degbar.data()
set.degbar.data <- function(geo2r_res) {
  deg_sum <- lapply(names(geo2r_res), function(geo2rname) {
    geo2r <- geo2r_res[[geo2rname]]
    deg_sum <- rep("b.Unperturbed", nrow(geo2r))
    deg_sum[which(geo2r[[grep('deg', colnames(geo2r))]] == -1)] <- "c.Downregulated"
    deg_sum[which(geo2r[[grep('deg', colnames(geo2r))]] ==  1)] <- "a.Upregulated"
    data.frame('dataset' = rep(geo2rname, nrow(geo2r)),
               'Regulation' = deg_sum)
  })
  Reduce(rbind, deg_sum)
}

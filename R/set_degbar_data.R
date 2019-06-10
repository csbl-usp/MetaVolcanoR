#' A function setting data format for DEG barplot visualization
#'
#' This function summarize the variable deg from the deg_def() function to
#' visualize as barplots the number of DEGs per input study
#' @param diffexp list of data.frame/data.table (s) output of the deg_def()
#'        function <list>
#' @keywords set DEG barplot data format
#' @return \code{data.frame} DEG by input
#' @export
#' @examples
#' data(diffexplist)
#' diffexp <- lapply(diffexplist, function(...) deg_def(..., "pvalue", 
#'            "Log2FC", 0.05, 0))
#' bardat <- set_degbar_data(diffexp)
#' head(bardat, 3)
set_degbar_data <- function(diffexp) {
  deg_sum <- lapply(names(diffexp), function(dename) {
    sdiffexp <- diffexp[[dename]]
    deg_sum <- rep("b.Unperturbed", nrow(sdiffexp))
    deg_sum[which(sdiffexp[[grep('deg', 
        colnames(sdiffexp))]] == -1)] <- "c.Downregulated"
    deg_sum[which(sdiffexp[[grep('deg', 
        colnames(sdiffexp))]] ==  1)] <- "a.Upregulated"
    data.frame('dataset' = rep(dename, nrow(sdiffexp)),
               'Regulation' = deg_sum)
  })
  Reduce(rbind, deg_sum)
}

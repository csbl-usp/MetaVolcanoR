#' A column renaming function for GEO2R outputs
#'
#' This function rename the probe and gene name columns from the GEO2R outputs
#' @param geo2r_list list of data.frame/data.table (s) with DE results where lines are genes
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param ncores the number of processors the user wants to use  <integer>
#' @keywords rename column
#' @export
#' @examples
#' rename.col()
rename.col <- function(geo2r_list, collaps, ncores) {
  mclapply(1:length(geo2r_list), function(nstudy) {
    geo2r <- geo2r_list[[nstudy]]
    colnames(geo2r) <- paste(colnames(geo2r), nstudy, sep = "_")
    if(collaps) {
      colnames(geo2r)[7] <- "Gene.symbol"
    } else {
      colnames(geo2r)[1] <- "Probe"
    }
    return(geo2r)
  }, mc.cores = ncores)
}	

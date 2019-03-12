#' A column renaming function for GEO2R outputs
#'
#' This function rename the probe and gene name columns from the GEO2R outputs
#' @param geo2r_list list of data.frame/data.table (s) with DE results where lines are genes
#' @param genenamecol the column name of the gene name variable <string>
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript variable <string>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param ncores the number of processors the user wants to use  <integer>
#' @keywords rename column
#' @export
#' @examples
#' rename.col()
rename.col <- function(geo2r_list, genenamecol, geneidcol, collaps, ncores) {
  ns <- names(geo2r_list)
  geo2rs <- mclapply(seq(geo2r_list), function(nstudy) {
    geo2r <- geo2r_list[[nstudy]]
    colnames(geo2r) <- paste(colnames(geo2r), nstudy, sep = "_")
    if(collaps) {
      colnames(geo2r)[grep(genenamecol, colnames(geo2r))] <- genenamecol
    } else {
      colnames(geo2r)[grep(geneidcol, colnames(geo2r))] <- geneidcol
    }
    return(geo2r)
  }, mc.cores = ncores)
  names(geo2rs) <- ns
  return(geo2rs)
}

#' A function to filter out probes standing for the same gene which have
#' contradictory DE direction  
#' 
#' This function check if the probes standing for the same gene has a unique 
#' value for the variable deg from deg.def() 
#' @param geo2r data.frame/data.table output of the deg.def() function
#' @param genenamecol the column name of the gene name variable <string>
#' @param pcriteria the column name of the Pval criteria to consider <string>
#' @keywords Probe collapse
#' @export
#' @examples
#' collapse.deg()
collapse.deg <- function(geo2r, genenamecol, pcriteria) {
  geo2r %>%
    arrange(!!as.name(pcriteria)) %>%
    filter(!duplicated(!!as.name(genenamecol)))
}

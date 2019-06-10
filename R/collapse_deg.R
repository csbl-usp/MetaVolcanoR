#' A function to filter out geneIDs standing for the same gene name  
#' 
#' This function to remove redundant geneIDs standing for the same gene name
#' @param diffexp data.frame/data.table output of the deg.def() function
#' @param genenamecol the column name of the gene name variable <string>
#' @param pcriteria the column name of the pvalue criteria to consider <string>
#' @keywords Probe collapse
#' @return \code{data.table} differential expression results with unique gene
#'         names
#' @export
#' @examples
#' data(diffexplist)
#' diffexp <- collapse_deg(diffexplist[[1]], "Symbol", "pvalue")
#' head(diffexp, 3)
collapse_deg <- function(diffexp, genenamecol, pcriteria) {
  diffexp %>%
    arrange(!!as.name(pcriteria)) %>%
    dplyr::filter(!duplicated(!!as.name(genenamecol)))
}

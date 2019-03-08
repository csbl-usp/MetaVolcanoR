#' A function to filter out probes standing for the same gene which have
#' contradictory DE direction  
#' 
#' This function check if the probes standing for the same gene has a unique 
#' value for the variable deg from deg.def() 
#' @param geo2r data.frame/data.table output of the deg.def() function
#' @keywords Probe collapse
#' @export
#' @examples
#' collapse.deg()
collapse.deg <- function(geo2r) {
  ugen <- names(which(table(geo2r[['Gene.symbol']]) == 1))
  dgen <- names(which(table(geo2r[['Gene.symbol']]) > 1))
  if(length(dgen) == 0) {
    geo2r
  } else {
    sdgen <- filter(geo2r, Gene.symbol %in% dgen)
    expdir <- summarize(group_by(sdgen, Gene.symbol), deg_sum = length(unique(deg)))
    rbind(filter(geo2r, Gene.symbol %in% ugen), 
          filter(sdgen, Gene.symbol %in% filter(expdir, deg_sum == 1)[['Gene.symbol']]))
  }
}

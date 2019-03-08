#' A function to draw the 'Vote-counting meta-analysis' MetaVolcano 
#'
#' This function draws the 'Vote-counting meta-analysis' MetaVolcano
#' @param meta_geo2r data.frame/data.table containing all the inputed GEO2R outputs
#' @param metathr the proportion of studies where a gene has to be consistently DE to be highlighted 'vote-counting meta-analysis'
#' @param nstud the number of inputed GEO2R outputs  <integer>
#' @param jobname name of the running job <string>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param outputfolder /path where to write the results/ <string>
#' @keywords write 'vote-counting meta-analysis' metavolcano
#' @export
#' @examples
#' draw.metavolcano()
draw.metavolcano <- function(meta_geo2r, metathr, nstud, jobname, collaps, outputfolder) {
  meta_geo2r <- mutate(meta_geo2r, 
                       dcol_vote = ifelse(ndeg >= nstud*metathr & ddeg >= nstud*metathr, "Up-regulated",
                                          ifelse(ndeg >= nstud*metathr & ddeg <= -nstud*metathr, "Down-regulated", 
                                                 "Unperturbed")))
  print(head(meta_geo2r,2))
  meta_geo2r[['Gene.symbol']] <- apply(select(meta_geo2r, matches('symbol')), 1, function(g) {
    paste(unique(g[!is.na(g)]), collapse = '///')
  })
  # --- Drawing volcano ggplotly 
  gg <- draw.mv.gplotly(meta_geo2r, nstud, metathr, collaps, FALSE)
  # --- Writing html device for offline visualization
  htmlwidgets::saveWidget(as_widget(gg), paste0(outputfolder, jobname, ".html"))
  meta_geo2r
}
#' A function to draw the 'Combining meta-analysis' MetaVolcano 
#'
#' This function draws the 'Combining meta-analysis' MetaVolcano
#' @param meta_geo2r data.frame/data.table containing all the inputed GEO2R outputs
#' @param pcriteria the column name of the Pval criteria to consider c("adj.P.Val", "P.Value") <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param metathr the proportion of top meta-DEGs to be highlighted 'combining meta-analysis'<double>
#' @param nstud the number of inputed GEO2R outputs  <integer>
#' @param jobname name of the running job <string>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param metap method for Pval combination/merge. c("Fisher")  <string>
#' @param metafc method for summarizing gene fold-changes across studies c("Mean", "Median") <string>
#' @param outputfolder /path where to write the results/ <string>
#' @keywords write 'combining meta-analysis' metavolcano
#' @export
#' @examples
#' draw.metavolcano.metap()
draw.metavolcano.metap <- function(meta_geo2r, pcriteria, genenamecol, foldchangecol, metathr, nstud, jobname, collaps, metap, metafc, outputfolder) {
  if(metap == "Fisher") {
    meta_geo2r <- mutate(meta_geo2r, 
                         metap = apply(select(meta_geo2r, matches(pcriteria)), 1, 
                                       function(p) { 
                                         pp <- as.numeric(p[which(!is.na(p))])
                                         if(length(pp) == 1) {
                                           pp
                                         } else {
                                           metap::sumlog(pp)$p
                                         }
                                       }))
  }
  if(metafc == "Mean") {
    meta_geo2r <- mutate(meta_geo2r,
                         metafc = apply(select(meta_geo2r, matches(foldchangecol)), 1, 
                                        function(...) mean(as.numeric(...), na.rm = TRUE)))
  } else if (metafc == "Median") {
    meta_geo2r <- mutate(meta_geo2r, 
                         metafc = apply(select(meta_geo2r, matches(foldchangecol)), 1, 
                                        function(...) median(as.numeric(...), na.rm = TRUE)))
  }
  
  meta_geo2r <- mutate(meta_geo2r, 
                       dcol_combin = ifelse(metap <= quantile(as.numeric(metap), 1-metathr, na.rm = TRUE) & 
                                              metafc >= quantile(as.numeric(metafc), metathr, na.rm = TRUE), "Up-regulated",
                                            ifelse(metap <= quantile(as.numeric(metap), 1-metathr, na.rm = TRUE) & 
                                                     metafc <= quantile(as.numeric(metafc), 1-metathr, na.rm = TRUE), "Down-regulated", 
                                                   "Unperturbed")))
  
  meta_geo2r[['Gene.symbol']] <- apply(select(meta_geo2r, matches(genenamecol)), 1, function(g) {
    paste(unique(g[!is.na(g)]), collapse = '///')
  })
  print(head(meta_geo2r))
  # --- Drawing volcano ggplotly 
  gg <- draw.mv.gplotly(meta_geo2r, nstud, metathr, collaps, TRUE)
  # --- Writing html device for offline visualization
  htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder), 
                                                '/combining_method_MetaVolcano_', jobname, ".html"))
  meta_geo2r
}

#' A function to run the 'Random Effect Model' MetaVolcano 
#'
#' This function runs the 'Random Effect Model' MetaVolcano section
#' @param geo2r_res list of data.frame/data.table (s) with DE results where lines are genes
#' @param pcriteria the column name of the Pval criteria to consider c("adj.P.Val", "P.Value") <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript variable <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param pvalue the Pval to use as threshold c(0:1) <double>
#' @param logfc the foldchange to use as DE threshold c(-Inf: Inf) <double>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param ncores the number of processors the user wants to use <integer>
#' @param cvar weather or not to calculate gene variance from confidence interval limits <logical>
#' @keywords write 'combining meta-analysis' metavolcano
#' @export
#' @examples
#' do.metafor()
do.metafor <- function(geo2r_res, pcriteria, genenamecol, geneidcol, foldchangecol, pvalue, logfc, collaps, jobname, outputfolder, ncores, cvar) {
  # --- Defining DEGs, criteria = Pmethod, pvalue and lofFC-value
  geo2r_res <- lapply(geo2r_res, function(...) deg.def(..., pcriteria, foldchangecol, pvalue, logfc))
  
  # --- Removing non-named genes
  #geo2r_res <- mclapply(geo2r_res, function(...) filter(..., Gene.symbol != ""), mc.cores = ncores)
  
  # ---- Calculating variance from coifidence intervals
  if(cvar == TRUE) {
    geo2r_res <- lapply(geo2r_res, function(...) calc.vi(...))
  }
  
  # --- Collapsing genes whose probes do not have the same expression pattern
  if(collaps) {
    geo2r_res_col <- lapply(geo2r_res, function(...) collapse.deg(..., genenamecol))
    meta_geo2r <- Reduce(function(x, y) merge(x, y, by = genenamecol, all = TRUE), geo2r_res)
  } else {
    meta_geo2r <- Reduce(function(x, y) merge(x, y, by = geneidcol, all = TRUE), geo2r_res)
  }
  
  #geo2r_res_col <- mclapply(geo2r_res_col, function(...) filter(..., !duplicated(Gene.symbol)), 
  #			  mc.cores = ncores)
  
  #names(geo2r_res_col) <- names(geo2r_res)
  
  #print(str(geo2r_res_col))
  
  meta_geo2r[['Gene.symbol']] <- apply(select(meta_geo2r, matches(genenamecol)), 1, function(g) {
    paste(unique(g[!is.na(g)]), collapse = '///')
  })
  
  
  meta_res <- cbind(meta_geo2r, Reduce(rbind, apply(meta_geo2r, 1, function(...) rnmodel(..., foldchangecol))))
  draw.metav(meta_res, jobname, outputfolder)
  meta_res
}

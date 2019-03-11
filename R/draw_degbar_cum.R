#' A function for the initial MetaVolcano overview
#'
#' This function runs the first Meta-Volcano section 
#' i) DEGs per study. Write .html barplot interactive plotly 
#' ii) Inverse-cummulative DEG distribution along the inputed GEO2R outputs. Write .html lineplot interactive plotly 
#' iii) Merge the inputed GEO2R outputs in one data.frame/data.table. Output the meta_geo2r data.frame object
#' @param geo2r_res list of data.frame/data.table (s) with DE results where lines are genes
#' @param pcriteria the column name of the Pval criteria to consider <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param pvalue the Pval to use as threshold c(0:1) <double>
#' @param logfc the foldchange to use as DE threshold c(-Inf: Inf) <double>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param ncores the number of processors the user wants to use <integer>
#' @keywords write metavolcano
#' @export
#' @examples
#' draw.degbar.cum()
draw.degbar.cum <- function(geo2r_res, pcriteria, foldchangecol, genenamecol, pvalue, logfc, collaps, jobname, outputfolder, ncores) {
  nstud <- length(geo2r_res)
  # --- Defining DEGs, criteria = Pmethod, pvalue and lofFC-value
  geo2r_res <- lapply(geo2r_res, function(...) deg.def(..., pcriteria, foldchangecol, pvalue, logfc))
  
  if (collaps) {
    # --- Removing non-named genes
    geo2r_res <- mclapply(geo2r_res, function(...) filter(..., !!as.name(genenamecol) != ""), mc.cores = ncores)
    
    # --- Collapsing genes whose probes do not have the same expression pattern
    geo2r_res_col <- lapply(geo2r_res, function(...) collapse.deg(..., genenamecol))
    #geo2r_res_col <- rename.col(geo2r_res_col, collaps, ncores)
    geo2r_res_col <- mclapply(geo2r_res_col, function(...) filter(..., !duplicated(!!as.name(genenamecol))), 
                              mc.cores = ncores)
    names(geo2r_res_col) <- names(geo2r_res)
    
    # --- Drawing DEGs by dataset
    gg <- draw.degbar(set.degbar.data(geo2r_res_col))
    # --- Writing html device for offline visualization
    htmlwidgets::saveWidget(as_widget(gg), paste0(outputfolder, "degbar_", jobname, ".html"))
    
    # --- merging DEG results	
    meta_geo2r <- Reduce(function(...) merge(..., by = genenamecol, all = TRUE), geo2r_res_col)
    
    # --- Defining new vars for visualization
    meta_geo2r[['ndeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum((r^2), na.rm = T))
    meta_geo2r[['ddeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum(r, na.rm = TRUE))
    
    # --- Drawing cDEGs by dataset
    gg <- draw.cum.freq(meta_geo2r, nstud)
    print(nstud)
    # --- Writing html device for offline visualization
    htmlwidgets::saveWidget(as_widget(gg), paste0(outputfolder, "cumdeg_", jobname, ".html"))
    
    #return(filter(meta_geo2r, ndeg != 0))
    return(meta_geo2r)
    
  } else {
    
    # --- Drawing DEGs by dataset
    gg <- draw.degbar(set.degbar.data(geo2r_res))
    # --- Writing html device for offline visualization
    htmlwidgets::saveWidget(as_widget(gg), paste0(outputfolder, "degbar_", jobname, ".html"))
    
    # --- merging DEG results	
    meta_geo2r <- Reduce(function(x, y) merge(x, y, by = 'Probe', all = TRUE), rename.col(geo2r_res, collaps, ncores))
    
    # --- Defining new vars for visualization
    meta_geo2r[['ndeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum((r^2), na.rm = T))
    meta_geo2r[['ddeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum(r, na.rm = TRUE))
    
    # --- Drawing cDEGs by dataset
    gg <- draw.cum.freq(meta_geo2r, nstud)
    # --- Writing html device for offline visualization
    htmlwidgets::saveWidget(as_widget(gg), paste0(outputfolder, "cumdeg_", jobname, ".html"))
    
    #		return(filter(meta_geo2r, ndeg != 0))
    return(meta_geo2r)
  }	
  
}

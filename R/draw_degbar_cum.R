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
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript variable <string>
#' @param pvalue the Pval to use as threshold c(0:1) <double>
#' @param logfc the foldchange to use as DE threshold c(-Inf: Inf) <double>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param draw wheather or not to draw the .html visualization 
#' @param ncores the number of processors the user wants to use <integer>
#' @keywords write metavolcano
#' @export
#' @examples
#' draw.degbar.cum()
draw.degbar.cum <- function(geo2r_res, pcriteria, foldchangecol, genenamecol, geneidcol, pvalue, logfc, collaps, jobname, outputfolder, draw, ncores) {
  nstud <- length(geo2r_res)
  # --- Defining DEGs
  geo2r_res <- lapply(geo2r_res, function(...) deg.def(..., pcriteria, foldchangecol, pvalue, logfc))
  
  if (collaps) {
    # --- Removing non-named genes
    geo2r_res <- mclapply(geo2r_res, function(g) {
      g %>%
        filter(!!as.name(genenamecol) != "") %>%
        filter(!is.na(!!as.name(genenamecol)))
      }, mc.cores = ncores)
    
    # --- Collapsing redundant geneIDs. Rataining the geneID with the smallest pcriteria
    geo2r_res <- mclapply(geo2r_res, function(g) {
      collapse.deg(g, genenamecol, pcriteria)
    }, mc.cores = ncores)
    
    # --- Drawing DEGs by dataset
    if(draw) {
      bardat <- set.degbar.data(geo2r_res)
      gg <- draw.degbar(bardat)
      # --- Writing html device for offline visualization
      htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder), "/deg_by_study_", jobname, ".html"))
    }
    
    # --- merging DEG results
    geo2r_res <- rename.col(geo2r_res, genenamecol, geneidcol, collaps, ncores)
    meta_geo2r <- Reduce(function(...) merge(..., by = genenamecol, all = TRUE), geo2r_res)
    
    # --- Defining new vars for visualization
    meta_geo2r[['ndeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum((r^2), na.rm = T))
    meta_geo2r[['ddeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum(r, na.rm = TRUE))
    
    # --- Drawing cDEGs by dataset
    if(draw) {
      gg <- draw.cum.freq(meta_geo2r, nstud)
      # --- Writing html device for offline visualization
      htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder), "/deg_InvCumDist_", jobname, ".html"))
    }
    
    # Return genes that were DE in at least one study
    return(filter(meta_geo2r, ndeg != 0))
  } else {
    gid <- sapply(geo2r_res, function(g) length(unique(g[[geneidcol]])) == nrow(g))
    if(all(gid)) {
      # --- Drawing DEGs by dataset
      if(draw) {
        bardat <- set.degbar.data(geo2r_res)
        gg <- draw.degbar(bardat)
        # --- Writing html device for offline visualization
        htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder), "/deg_by_study_", jobname, ".html"))
      }
      
      # --- merging DEG results	
      geo2r_res <- rename.col(geo2r_res, genenamecol, geneidcol, collaps, ncores)
      meta_geo2r <- Reduce(function(...) merge(..., by = geneidcol, all = TRUE), geo2r_res)
      
      # --- Defining new vars for visualization
      meta_geo2r[['ndeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum((r^2), na.rm = T))
      meta_geo2r[['ddeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum(r, na.rm = TRUE))
      
      # --- Drawing cDEGs by dataset
      if(draw) {
        gg <- draw.cum.freq(meta_geo2r, nstud)
        # --- Writing html device for offline visualization
        htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder), "/deg_InvCumDist_", jobname, ".html"))
      }
      return(filter(meta_geo2r, ndeg != 0))
    } else {
      stop("the geneidcol contains duplicated values, consider to set collaps=TRUE")
    }
  }
}

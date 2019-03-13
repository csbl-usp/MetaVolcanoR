#' A function to draw the 'Vote-counting meta-analysis' MetaVolcano 
#'
#' This function draws the 'Vote-counting meta-analysis' MetaVolcano
#' @param geo2r_res list of data.frame/data.table (s) with DE results where lines are genes
#' @param pcriteria the column name of the Pval criteria to consider <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript variable <string>
#' @param pvalue the Pval to use as threshold c(0:1) <double>
#' @param logfc the foldchange to use as DE threshold c(-Inf: Inf) <double>
#' @param metathr the proportion of studies a gene has to be DEG to be considered cDEG <double>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param draw wheather or not to draw the .html visualization 
#' @param ncores the number of processors the user wants to use <integer>
#' @keywords write 'vote-counting meta-analysis' metavolcano
#' @export
#' @examples
#' draw.metavolcano()
draw.metavolcano <- function(geo2r_res, pcriteria, foldchangecol, genenamecol, geneidcol, pvalue, logfc, metathr, collaps, jobname, outputfolder, draw, ncores) {
  nstud <- length(geo2r_res)
  # --- Subsetting DE inputs
  geo2r_res <- lapply(geo2r_res, function(...) select(..., matches(paste(c(pcriteria, foldchangecol, genenamecol, geneidcol), collapse = '|'))))
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
    
    
    # --- merging DEG results
    geo2r_res <- rename.col(geo2r_res, genenamecol, geneidcol, collaps, ncores)
    meta_geo2r <- Reduce(function(...) merge(..., by = genenamecol, all = TRUE), geo2r_res)
    
    # --- Defining new vars for visualization
    meta_geo2r[['ndeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum((r^2), na.rm = T))
    meta_geo2r[['ddeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum(r, na.rm = TRUE))
    
    # --- Defining perturbation based on the number of studies were a gene were identified as DE
    meta_geo2r <- mutate(meta_geo2r, 
                         dcol_vote = ifelse(ndeg >= nstud*metathr & ddeg >= nstud*metathr, "Up-regulated",
                                            ifelse(ndeg >= nstud*metathr & ddeg <= -nstud*metathr, "Down-regulated", 
                                                   "Unperturbed")))
    # --- 
    meta_geo2r[['Gene.symbol']] <- apply(select(meta_geo2r, matches(genenamecol)), 1, function(g) {
      paste(unique(g[!is.na(g)]), collapse = '///')
    })
    
    # --- Drawing cDEGs by dataset
    if(draw) {
      # --- Drawing volcano ggplotly 
      gg <- draw.mv.gplotly(meta_geo2r, nstud, metathr, collaps, FALSE, genenamecol)
      # --- Writing html device for offline visualization
      htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder), '/votecounting_metavolcano_', jobname, ".html"))
    }
    
    # Return genes that were highlighted as cDEG
    return(filter(meta_geo2r, dcol_vote != "Unperturbed"))
    
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
      
      # --- Defining perturbation based on the number of studies were a gene were identified as DE
      meta_geo2r <- mutate(meta_geo2r, 
                           dcol_vote = ifelse(ndeg >= nstud*metathr & ddeg >= nstud*metathr, "Up-regulated",
                                              ifelse(ndeg >= nstud*metathr & ddeg <= -nstud*metathr, "Down-regulated", 
                                                     "Unperturbed")))
      # --- 
      meta_geo2r[['Gene.symbol']] <- apply(select(meta_geo2r, matches(genenamecol)), 1, function(g) {
        paste(unique(g[!is.na(g)]), collapse = '///')
      })
      
      # --- Drawing cDEGs by dataset
      if(draw) {
        
        # --- Drawing volcano ggplotly 
        gg <- draw.mv.gplotly(meta_geo2r, nstud, metathr, collaps, FALSE,genenamecol)
        # --- Writing html device for offline visualization
        htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder), '/votecounting_metavolcano_', jobname, ".html"))
        
      }
      
      # Return genes that were highlighted as cDEG
      return(filter(meta_geo2r, dcol_vote != "Unperturbed"))
      
    } else {
      
      stop("the geneidcol contains duplicated values, consider to set collaps=TRUE")
      
    }
    
  }
  
}

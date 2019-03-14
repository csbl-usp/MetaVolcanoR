#' A function to draw the 'Combining meta-analysis' MetaVolcano
#'
#' This function draws the 'Combining meta-analysis' MetaVolcano
#' @param geo2r_res list of data.frame/data.table (s) with DE results where lines are genes
#' @param pcriteria the column name of the Pval criteria to consider c("adj.P.Val", "P.Value") <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript variable <string>
#' @param pvalue the Pval to use as threshold c(0:1) <double>
#' @param logfc the foldchange to use as DE threshold c(-Inf: Inf) <double>
#' @param metap method for Pval combination/merge. c("Fisher")  <string>
#' @param metafc method for summarizing gene fold-changes across studies c("Mean", "Median") <string>
#' @param metathr the proportion of top meta-DEGs to be highlighted 'combining meta-analysis'<double>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param draw wheather or not to draw the .html visualization <logical>
#' @param ncores the number of processors the user wants to use <integer>
#' @keywords write 'combining meta-analysis' metavolcano
#' @export
#' @examples
#' draw.metavolcano.metap()
draw.metavolcano.metap <- function(geo2r_res, pcriteria, foldchangecol, genenamecol, geneidcol, pvalue, logfc, metap, metafc, metathr, collaps, jobname, outputfolder, draw, ncores) {
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

    # --- Combining Pvalues with Fisher method
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

    # --- Combining FC values by either mean or median summary methods
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


    # --- Drawing cDEGs by dataset
    if(draw) {

      # --- Drawing volcano ggplotly
      gg <- draw.mv.gplotly(meta_geo2r, nstud, metathr, collaps, TRUE, genenamecol)
      # --- Writing html device for offline visualization
      htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder),
                                                    '/combining_method_MetaVolcano_', jobname, ".html"))

    }

    # Return genes that were highlighted as cDEG
    return(filter(meta_geo2r, dcol_combin != "Unperturbed"))

  } else {

    gid <- sapply(geo2r_res, function(g) length(unique(g[[geneidcol]])) == nrow(g))

    if(all(gid)) {

      # --- merging DEG results
      geo2r_res <- rename.col(geo2r_res, genenamecol, geneidcol, collaps, ncores)
      meta_geo2r <- Reduce(function(...) merge(..., by = geneidcol, all = TRUE), geo2r_res)

      # --- Combining Pvalues with Fisher method
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

      # --- Combining FC values by either mean or median summary methods
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

      # --- Drawing cDEGs by dataset
      if(draw) {

        # --- Drawing volcano ggplotly
        gg <- draw.mv.gplotly(meta_geo2r, nstud, metathr, collaps, TRUE, geneidcol)
        # --- Writing html device for offline visualization
        htmlwidgets::saveWidget(as_widget(gg), paste0(normalizePath(outputfolder),
                                                      '/combining_method_MetaVolcano_', jobname, ".html"))

      }

      # Return genes that were highlighted as cDEG
      return(filter(meta_geo2r, dcol_combin != "Unperturbed"))

    } else {

      stop("the geneidcol contains duplicated values, consider to set collaps=TRUE")

    }

  }

}

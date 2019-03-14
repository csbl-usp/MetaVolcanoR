#' A function to run the 'Random Effect Model' MetaVolcano
#'
#' This function runs the 'Random Effect Model' MetaVolcano section
#' @param geo2r_res list of data.frame/data.table (s) with DE results where lines are genes
#' @param pcriteria the column name of the Pval criteria to consider c("adj.P.Val", "P.Value") <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript variable <string>
#' @param pvalue the Pval to use as threshold c(0:1) <double>
#' @param logfc the foldchange to use as DE threshold c(-Inf: Inf) <double>
#' @param collaps if probes should be collapsed based on the DE direction <logical>
#' @param llcol left limit of the fold change coinfidence interval variable name <string>
#' @param rlcol right limit of the fold change coinfidence interval variable name <string>
#' @param vcol name of the fold change variance variable <string>
#' @param cvar weather or not to calculate gene variance from confidence interval limits <logical>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param draw wheather or not to draw the .html visualization <logical>
#' @param ncores the number of processors the user wants to use <integer>
#' @keywords write 'combining meta-analysis' metavolcano
#' @export
#' @examples
#' do.metafor()
do.metafor <- function(geo2r_res, pcriteria, foldchangecol, genenamecol, geneidcol, pvalue, logfc, collaps, llcol, rlcol, vcol, cvar, jobname, outputfolder, draw, ncores) {
  nstud <- length(geo2r_res)
  # --- Subsetting DE inputs
  geo2r_res <- lapply(geo2r_res, function(...) select(..., matches(paste(c(pcriteria, foldchangecol, genenamecol, geneidcol, llcol, rlcol, vcol), collapse = '|'))))
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

    # ---- Calculating variance from coifidence interval
    if(cvar == TRUE) {
      geo2r_res <- lapply(geo2r_res, function(...) calc.vi(..., llcol, rlcol))
    }

    # --- merging DEG results
    geo2r_res <- rename.col(geo2r_res, genenamecol, geneidcol, collaps, ncores)
    meta_geo2r <- Reduce(function(...) merge(..., by = genenamecol, all = TRUE), geo2r_res)

    # --- Defining new vars for visualization
    meta_geo2r[['ndeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum((r^2), na.rm = T))
    meta_geo2r[['ddeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum(r, na.rm = TRUE))

    meta_geo2r <- cbind(meta_geo2r, Reduce(rbind, apply(meta_geo2r, 1, function(...) rnmodel(..., foldchangecol))))

    # --- Drawing metavolcano
    if(draw) {
      draw.metav(meta_geo2r, jobname, outputfolder, genenamecol)
    }

    # Return REM results for all genes
    return(meta_geo2r)

  } else {

    gid <- sapply(geo2r_res, function(g) length(unique(g[[geneidcol]])) == nrow(g))

    if(all(gid)) {

      # ---- Calculating variance from coifidence interval
      if(cvar == TRUE) {
        geo2r_res <- lapply(geo2r_res, function(...) calc.vi(..., llcol, rlcol))
      }

      # --- merging DEG results
      geo2r_res <- rename.col(geo2r_res, genenamecol, geneidcol, collaps, ncores)
      meta_geo2r <- Reduce(function(...) merge(..., by = geneidcol, all = TRUE), geo2r_res)

      # --- Defining new vars for visualization
      meta_geo2r[['ndeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum((r^2), na.rm = T))
      meta_geo2r[['ddeg']] <- apply(select(meta_geo2r, matches("deg_")), 1, function(r) sum(r, na.rm = TRUE))

      meta_geo2r <- cbind(meta_geo2r, Reduce(rbind, apply(meta_geo2r, 1, function(...) rnmodel(..., foldchangecol))))

      # --- Drawing metavolcano
      if(draw) {
        draw.metav(meta_geo2r, jobname, outputfolder, geneidcol)
      }

      # Return REM results for all genes
      return(meta_geo2r)

    } else {

      stop("the geneidcol contains duplicated values, consider to set collaps=TRUE")

    }

  }

}

#' A function to draw the 'Combining meta-analysis' MetaVolcano
#'
#' This function draws the 'Combining meta-analysis' MetaVolcano
#' @param diffexp list of data.frame/data.table (s) with DE results where lines
#'        are genes
#' @param pcriteria the column name of the Pval criteria to consider 
#'        c("adj.P.Val", "P.Value") <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript
#'        variable <string>
#' @param metafc method for summarizing gene fold-changes across studies 
#'        c("Mean", "Median") <string>
#' @param metathr top percentage of perturbed genes to be highlighted <double>
#' @param collaps if probes should be collapsed based on the DE direction
#'        <logical>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param draw wheather or not to draw the .pdf or .html visualization 
#'        <c(NULL, "PDF", "HTML")>
#' @param ncores the number of processors the user wants to use <integer>
#' @keywords write 'combining meta-analysis' metavolcano
#' @export
#' @examples
#' combining_mv()
combining_mv <- function(diffexp=list(), pcriteria="pvalue", 
			 foldchangecol="Log2FC", genenamecol="Symbol", 
			 geneidcol=NULL, metafc="Mean", metathr=0.01, 
			 collaps="FALSE", jobname="MetaVolcano", 
			 outputfolder=".", draw="HTML", ncores=1) {
    	
    if (collaps) {
        # --- Removing non-named genes
        diffexp <- mclapply(diffexp, function(g) {
            g %>%
            dplyr::filter(!!as.name(genenamecol) != "") %>%
            dplyr::filter(!is.na(!!as.name(genenamecol)))
        }, mc.cores = ncores)

        # --- Collapsing redundant geneIDs. Rataining the geneID with the 
        # --- smallest pcriteria
        diffexp <- mclapply(diffexp, function(g) {
            collapse_deg(g, genenamecol, pcriteria)
        }, mc.cores = ncores)

	# --- Subsetting the diffexp inputs
	diffexp <- lapply(diffexp, function(...) dplyr::select(...,
			  dplyr::matches(paste(c(pcriteria, foldchangecol,
				genenamecol), collapse = '|'))))

        # --- merging DEG results
        diffexp <- rename_col(diffexp, genenamecol, ncores)
        meta_diffexp <- Reduce(function(...) merge(..., by = genenamecol, 
						   all = TRUE), diffexp)
	genecol <- genenamecol

    } else {

	if(is.null(geneidcol)) {
	    geneidcol <- genenamecol
	}

	# Testing if geneIDs are unique
	gid <- vapply(diffexp, function(g) {
                    length(unique(g[[geneidcol]])) == nrow(g)
                }, 
	    logical(1))
        
        if(all(gid)) {
        
	    # --- Subsetting the diffexp inputs
            diffexp <- lapply(diffexp, function(...) dplyr::select(...,
			  dplyr::matches(paste(c(pcriteria, foldchangecol,
				geneidcol), collapse = '|'))))

            # --- merging DEG results	
            diffexp <- rename_col(diffexp, geneidcol, ncores)
            meta_diffexp <- Reduce(function(...) merge(..., 
						       by = geneidcol, 
						       all = TRUE), diffexp)

	    genecol <- geneidcol

	} else {
		
	    stop("the geneidcol contains duplicated values, consider to 
		 set collaps=TRUE")

	}
    }    	
    
    # --- Combining Pvalues with Fisher method
    meta_diffexp <- dplyr::mutate(meta_diffexp,
        metap = apply(dplyr::select(meta_diffexp, dplyr::matches(pcriteria)), 1,
            function(p) {
	        pp <- as.numeric(p[which(!is.na(p))])
	        if(length(pp) == 1) {
		    pp
		} else {
		    metap::sumlog(pp)$p
		}
	}))
    
    # --- Combining fold-change by either mean or median summary methods
    if(metafc == "Mean") {
        meta_diffexp <- dplyr::mutate(meta_diffexp,
	    metafc = apply(dplyr::select(meta_diffexp,
					 dplyr::matches(foldchangecol)), 1,
                                          function(...) mean(as.numeric(...), 
							     na.rm = TRUE)))
    } else if (metafc == "Median") {
        meta_diffexp <- dplyr::mutate(meta_diffexp,
            metafc = apply(dplyr::select(meta_diffexp,
					 dplyr::matches(foldchangecol)), 1,
                                          function(...) median(as.numeric(...),
							       na.rm = TRUE)))
    } else {
	    stop("Oops! Please check the provided metafc parameter. Try either
		 Mean or Median")
    }
    

    # Highlighting the top perturbed genes
    meta_diffexp <- meta_diffexp %>%
        dplyr::mutate(idx = metafc*-log10(metap))
    meta_diffexp <- meta_diffexp %>%
	dplyr::mutate(degcomb = ifelse(idx < quantile(meta_diffexp[['idx']], 
						      metathr/2), 
				       '0.Down-regulated',
				ifelse(idx > quantile(meta_diffexp[['idx']], 
						      1-(metathr/2)), 
				       '2.Up-regulated', '1.Unperturbed'))) %>%
    dplyr::arrange(-abs(idx))

    # --- Drawing combining MetaVolcano
    gg <- plot_mv(meta_diffexp, NULL, genecol, TRUE, metafc)
    
    if(draw == "HTML") {

        # --- Writing html device for offline visualization
        htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	        '/combining_method_MetaVolcano_', jobname, ".html"))

    } else if(draw == "PDF") {

        # --- Writing PDF visualization
        pdf(paste0(normalizePath(outputfolder), 
            '/combining_method_MetaVolcano_', jobname, ".pdf"), 
	     width = 4, height = 5)
	        plot(gg)
        dev.off()

    } else {
		
        stop("Seems like you did not provide a right 'draw' parameter. 
              Try NULL, 'PDF' or 'HTML'")

    }
    
    # Set combining result
    icols <- paste(c(genecol, pcriteria, foldchangecol), collapse="|")
    rcols <- paste(c(genecol, "metafc", "metap", "idx"), collapse="|")
    result <- new('MetaVolcano', 
		  input=dplyr::select(meta_diffexp, 
				      dplyr::matches(icols)),
		  inputnames=names(diffexp),
		  metaresult=dplyr::select(meta_diffexp,
				       dplyr::matches(rcols)),
		  MetaVolcano=gg,
		  degfreq=ggplot()
		  )
    return(result)	
}

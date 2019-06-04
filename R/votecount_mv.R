#' A function to draw the 'Vote-counting meta-analysis' MetaVolcano
#'
#' This function draws the vote-counting meta-analysis MetaVolcano
#' @param diffexp list of data.frame/data.table (s) with DE results where lines
#'        are genes
#' @param pcriteria the column name of the Pval criteria to consider <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript 
#'        variable <string>
#' @param pvalue the Pval to use as threshold c(0:1) <double>
#' @param foldchange the foldchange to use as DE threshold c(-Inf: Inf) <double>
#' @param metathr the proportion of studies a gene has to be DEG to be 
#'        considered cDEG <double>
#' @param collaps if probes should be collapsed based on the DE 
#'        direction <logical>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param draw wheather or not to draw a .pdf or .html visualization 
#'        <c(NULL, 'PDF', 'HTML')>
#' @param ncores the number of processors the user wants to use <integer>
#' @keywords write 'vote-counting meta-analysis' metavolcano
#' @export
#' @examples
#' votecount_mv()
votecount_mv <- function(diffexp=list(), pcriteria="pvalue", 
			      foldchangecol="Log2FC", genenamecol="Symbol", 
			      geneidcol=NULL, pvalue=0.05, foldchange=0, 
			      metathr=0.01, collaps=FALSE, 
			      jobname="MetaVolcano", outputfolder=".", 
			      draw="HTML", ncores=1) {

    nstud <- length(diffexp)
  
    # --- Defining DEGs
    diffexp <- lapply(diffexp, function(...) deg_def(..., pcriteria, 
						     foldchangecol, 
						     pvalue, foldchange))

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
				genenamecol, '^deg$'), collapse = '|'))))

	# Setting data for DEG by study visualization
	bardat <- set_degbar_data(diffexp)

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
				geneidcol, '^deg$'), collapse = '|'))))
            # DEG by study data setting
	    bardat <- set_degbar_data(diffexp)

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
    # --- Defining new vars for visualization
    meta_diffexp[['ndeg']] <- apply(dplyr::select(meta_diffexp, 
					   dplyr::matches("deg_")), 
			       1, function(r) sum((r^2), na.rm = TRUE))
    meta_diffexp[['ddeg']] <- apply(dplyr::select(meta_diffexp, 
					   dplyr::matches("deg_")), 
			       1, function(r) sum(r, na.rm = TRUE))

    # Highlighting the top perturbed genes
    meta_diffexp <- meta_diffexp %>%
        dplyr::mutate(idx = ddeg*ndeg)
    meta_diffexp <- meta_diffexp %>%
        dplyr::mutate(degvcount = ifelse(idx < quantile(meta_diffexp[['idx']], 
						      metathr/2), 
				       '0.Down-regulated',
				ifelse(idx > quantile(meta_diffexp[['idx']], 
						      1-(metathr/2)), 
				       '2.Up-regulated', '1.Unperturbed')))

    # --- Drawing DEGs by dataset
    if(!is.null(draw)) {

	gg <- draw_degbar(bardat)
	ff <- draw_cum_freq(meta_diffexp, nstud)
        mv <- plot_mv(meta_diffexp, nstud, genecol, FALSE, NULL)

	if(draw == "HTML") {
        
	    # --- Writing html device for offline visualization
	    htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
		paste0(normalizePath(outputfolder), "/deg_by_study_", 
		       jobname, ".html"))
	    htmlwidgets::saveWidget(as_widget(ggplotly(ff)), 
		paste0(normalizePath(outputfolder), "/deg_InvCumDist_", 
		       jobname, ".html"))
	    htmlwidgets::saveWidget(as_widget(ggplotly(mv)), 
		paste0(normalizePath(outputfolder), 
		       '/votecounting_metavolcano_', jobname, ".html"))

	} else if(draw == "PDF") {

	    # --- Writing PDF visualization
	    pdf(paste0(normalizePath(outputfolder),
		       "/deg_by_study_", jobname,
		       ".pdf"), width = 7, height = 4)
	        plot(gg)
	    dev.off()

	    pdf(paste0(normalizePath(outputfolder),
		       "/deg_InvCumDist_", jobname,
		       ".pdf"), width = 7, height = 4)
	        plot(ff)
	    dev.off()

	    pdf(paste0(normalizePath(outputfolder),
		       "/votecounting_metavolcano_", jobname,
		       ".pdf"), width = 7, height = 4)
	        plot(mv)
	    dev.off()

		
	} else {
		
	    stop("Seems like you did not provide a right 'draw' parameter. 
		 Try NULL, 'PDF' or 'HTML'")

    	}
    }
       
    # Return genes that were DE in at least one study
    return(dplyr::filter(meta_diffexp, ndeg != 0))
}

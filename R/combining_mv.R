#' @importFrom stats median quantile
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @importFrom plotly as_widget
#' @importFrom htmlwidgets saveWidget
NULL

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
#' @param render A boolean parameter that determines whether the plot should be rendered. 
#' If `TRUE`, the function will produce and save the plot based on the specified `draw` 
#' parameter (either as an HTML or PDF file). If `FALSE` (default), no plot will be 
#' rendered or saved. It's useful for cases where you might want to run the function 
#' for its side effects or calculations without necessarily visualizing the result.
#' @return \code{MetaVolcano} object
#' @keywords write 'combining meta-analysis' metavolcano
#' @export
#' @examples
#' data(diffexplist)
#' mv <- combining_mv(diffexplist)
#' str(mv)
combining_mv <- function(diffexp=list(), pcriteria="pvalue", 
			 foldchangecol="Log2FC", genenamecol="Symbol", 
			 geneidcol=NULL, metafc="Mean", metathr=0.01, 
			 collaps="FALSE", jobname="MetaVolcano", 
			 outputfolder=".", draw="HTML", render = F) {
    	
    if(!draw %in% c('PDF', 'HTML')) {
		
        stop("Oops! Seems like you did not provide a right 'draw' parameter. 
              Try 'PDF' or 'HTML'")

    } else if(!metafc %in% c('Mean', 'Median')) {

	    stop("Oops! Please check the provided metafc parameter. Try either
		 Mean or Median")
    }

    if (collaps) {
        # --- Removing non-named genes
        diffexp <- lapply(diffexp, function(g) {
            g %>%
            dplyr::filter(!!as.name(genenamecol) != "") %>%
            dplyr::filter(!is.na(!!as.name(genenamecol))) %>%
            dplyr::filter(!!as.name(genenamecol) != "NA")
        })

        # --- Collapsing redundant geneIDs. Rataining the geneID with the 
        # --- smallest pcriteria
        diffexp <- lapply(diffexp, function(g) {
            collapse_deg(g, genenamecol, pcriteria)
        })

	# --- Subsetting the diffexp inputs
	diffexp <- lapply(diffexp, function(...) dplyr::select(...,
			  dplyr::matches(paste(c(pcriteria, foldchangecol,
				genenamecol), collapse = '|'))))

        # --- merging DEG results
        diffexp <- rename_col(diffexp, genenamecol)
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
            diffexp <- rename_col(diffexp, geneidcol)
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
    dplyr::select(meta_diffexp, 
        dplyr::matches(paste(c(genecol, pcriteria), collapse = "|"))) %>%
	    tidyr::gather(study, value, -!!as.name(genecol)) %>%
	    dplyr::filter(!is.na(value)) %>%
	    dplyr::group_by(!!as.name(genecol)) %>%
	    dplyr::summarize('metap' = tryCatch({ 
	        metap::sumlog(value)$p},
		    error = function(e){ return(NA) })) %>%
	    dplyr::filter(!is.na(metap)) -> meta_p

    
    meta_diffexp <- merge(meta_diffexp, meta_p, by = genecol)    

    # --- Combining fold-change by either mean or median summary methods
    dplyr::select(meta_diffexp, 
        dplyr::matches(paste(c(genecol, foldchangecol), collapse = "|"))) %>%
	    tidyr::gather(study, value, -!!as.name(genecol)) %>%
	    dplyr::filter(!is.na(value)) %>%
	    dplyr::group_by(!!as.name(genecol)) -> meta_fc

    if(metafc == "Mean") {
	meta_fc %>%
	    dplyr::summarize('metafc' = mean(as.numeric(value), 
				      na.rm = TRUE)) -> meta_fc

    } else if (metafc == "Median") {

	meta_fc %>%
	    dplyr::summarize('metafc' = median(as.numeric(value), 
				      na.rm = TRUE)) -> meta_fc
    }

    meta_diffexp <- merge(meta_diffexp, meta_fc, by = genecol)
	
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
    if(render) {
      if(draw == "HTML") {
  
          # --- Writing html device for offline visualization
          saveWidget(as_widget(ggplotly(gg)), 
              paste0(normalizePath(outputfolder), 
  	        '/combining_method_MetaVolcano_', jobname, ".html"))
  
      } else if(draw == "PDF") {
  
          # --- Writing PDF visualization
          pdf(paste0(normalizePath(outputfolder), 
              '/combining_method_MetaVolcano_', jobname, ".pdf"), 
  	     width = 4, height = 5)
  	        plot(gg)
          dev.off()
  
      } 
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

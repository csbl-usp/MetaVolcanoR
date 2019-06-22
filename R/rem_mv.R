#' @importFrom parallel mclapply
#' @importFrom topconfects normal_confects
#' @importFrom methods new 'slot<-' show
#' @importFrom plotly as_widget ggplotly
#' @importFrom htmlwidgets saveWidget 
#' @import dplyr

setOldClass('gg')
setOldClass('ggplot')

#' An S4 class to represent MetaVolcanoR results
#' 
#' @slot input merged diiferential expression inputs \code{data.frame} 
#' @slot inputnames names of the differential expression inputs \code{character}
#' @slot metaresult meta-analysis results \code{data.frame}
#' @slot MetaVolcano plot with meta-analysis results
#' @slot degfreq supplementary figure of the vote-counting MetaVolcano

setClass('MetaVolcano', slots = list(input='data.frame',
				     inputnames='character',
				     metaresult='data.frame',
				     MetaVolcano='gg',
				     degfreq='gg'
				     ))

#' A function to perform the Random Effect Model (REM) MetaVolcano
#'
#' This function runs the 'Random Effect Model' MetaVolcano section
#' @param diffexp list of data.frame/data.table (s) with DE results where lines
#'        are genes
#' @param pcriteria the column name of the pvalue variable <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param genenamecol the column name of the gene name variable <string>
#' @param geneidcol the column name of the gene ID/probe/oligo/transcript
#'        variable <string>
#' @param collaps if probes should be collapsed based on the DE direction
#'        <logical>
#' @param llcol left limit of the fold change coinfidence interval variable
#'        name <string>
#' @param rlcol right limit of the fold change coinfidence interval variable
#'        name <string>
#' @param vcol name of the fold change variance variable <string>
#' @param cvar weather or not to calculate gene variance from confidence 
#'        interval limits <logical>
#' @param metathr top percentage of perturbed genes to be highlighted <double>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/
#' @param draw wheather or not to draw the .html visualization <logical>
#' @param ncores the number of processors the user wants to use <integer>
#' @keywords write 'combining meta-analysis' metavolcano
#' @return MetaVolcano object
#' @export
#' @examples
#' data(diffexplist)
#' diffexplist <- lapply(diffexplist, function(del) {
#'     dplyr::filter(del, grepl("MP", Symbol))
#' })
#' mv <- rem_mv(diffexplist, metathr = 0.1)
#' str(mv)
rem_mv <- function(diffexp=list(), pcriteria="pvalue", foldchangecol="Log2FC",
		   genenamecol="Symbol", geneidcol=NULL, collaps=FALSE, 
		   llcol="CI.L", rlcol="CI.R", vcol=NULL, cvar=TRUE, 
		   metathr=0.01, jobname="MetaVolcano", outputfolder=".", 
		   draw='HTML', ncores=1) {

    
    if(!draw %in% c('PDF', 'HTML')) {
		
        stop("Oops! Seems like you did not provide a right 'draw' parameter. 
              Try 'PDF' or 'HTML'")

    }
    
    # ---- Calculating variance from coifidence interval
    if(cvar == TRUE) {
        diffexp <- lapply(diffexp, function(...) calc_vi(..., llcol, rlcol))
    	vcol <- 'vi'
    } else {
    	if(is.null(vcol)) {
	    stop("Oops! If cvar=FALSE, you should provide a variance stimate, 
		  Please, check the vcol parameter.")
	}
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
			  dplyr::matches(paste(c(genenamecol, foldchangecol,
			                         llcol, rlcol, vcol), 
					       collapse = '|'))))

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
			  dplyr::matches(paste(c(geneidcol, foldchangecol,
						 llcol, rlcol, vcol),
					       collapse = '|'))))

            # --- merging the diffexp inputs	
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
	
    # Calculating the REM summary (metafor)
    # computational intensive parallel run recommended
    meta_diffexp <- cbind(meta_diffexp,
                        do.call(rbind,
                                mclapply(split(meta_diffexp, 
					      meta_diffexp[[genecol]]),
                                        function(g) {
                                            remodel(g, foldchangecol, vcol)
                                        }, mc.cores = ncores)
                        )
    )

    # --- Topconfects ranking
    meta_diffexp <- meta_diffexp %>%
        dplyr::mutate(se = (randomCi.ub - randomCi.lb)/3.92) %>% # 95% conf.int
        dplyr::mutate(index = seq(nrow(meta_diffexp)))

    confects <- normal_confects(meta_diffexp$randomSummary, 
				se=meta_diffexp$se, 
				fdr=0.05, 
				full=TRUE)

    meta_diffexp <- merge(meta_diffexp, 
			  dplyr::select(confects$table, c(index, `rank`)), 
			  by = 'index', all = TRUE)
    
    meta_diffexp <- dplyr::arrange(meta_diffexp, `rank`)
   
    # --- Draw REM MetaVolcano
    gg <- plot_rem(meta_diffexp, jobname, outputfolder, genecol, metathr)

    if(draw == "HTML") {
        
        # --- Writing html device for offline visualization
        saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/RandomEffectModel_MetaVolcano_", 
	           jobname, ".html"))

    } else if(draw == "PDF") {

        # --- Writing PDF visualization
	pdf(paste0(normalizePath(outputfolder),
	           "/RandomEffectModel_MetaVolcano_", jobname,
	           ".pdf"), width = 7, height = 6)
	     plot(gg)
	dev.off()

    } 

    # Set REM result
    icols <- paste(c(genecol, pcriteria, foldchangecol, llcol, rlcol, vcol), 
		   collapse="|^")
    rcols <- paste(c(genecol, "^random", "^het_", "^error$", "^rank$", 
                   "signcon"), collapse="|")
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

#' @importFrom parallel mclapply
#' @importFrom topconfects normal_confects
#' @importFrom methods new 'slot<-' show
#' @importFrom plotly as_widget ggplotly
#' @importFrom htmlwidgets saveWidget
#' @import dplyr
NULL

#' Random Effect Model (REM) MetaVolcano
#'
#' This function performs Random Effect Model meta-analysis on differential expression results.
#' This function runs the 'Random Effect Model' MetaVolcano section
#' @export
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
#' @param colors named vector of colors for gradient: c(low, mid, high, na)
#' @param point_size size of points in the plot
#' @param label_genes character vector of specific genes to label
#' @param label_top_n number of top genes (by p-value) to label automatically
#' @param label_size size of gene name labels
#' @param plot_title custom plot title (NULL for no title)
#' @param show_legend whether to show the legend (TRUE/FALSE)
#' #' @param render A boolean parameter that determines whether the plot should be rendered. 
#' If `TRUE`, the function will produce and save the plot based on the specified `draw` 
#' parameter (either as an HTML or PDF file). If `FALSE` (default), no plot will be 
#' rendered or saved. It's useful for cases where you might want to run the function 
#' for its side effects or calculations without necessarily visualizing the result.
#' @keywords write 'combining meta-analysis' metavolcano
#' @return MetaVolcano object
#' @export
#' @examples
#' data(diffexplist)
#' diffexplist <- lapply(diffexplist, function(del) {
#'     dplyr::filter(del, grepl("MP", Symbol))
#' })
rem_mv <- function(diffexp=list(), pcriteria="pvalue", foldchangecol="Log2FC",
		   genenamecol="Symbol", geneidcol=NULL, collaps=FALSE,
		   llcol="CI.L", rlcol="CI.R", vcol=NULL, cvar=TRUE,
		   metathr=0.01, jobname="MetaVolcano", outputfolder = tempdir(),
		   draw='HTML', ncores=1,
		   colors = c(low = "#083e46", mid = "white", high = "#811820", na = "grey80"),
		   point_size = 0.6,
		   label_genes = NULL,
		   label_top_n = NULL,
		   label_size = 3,
		   plot_title = NULL,
		   show_legend = TRUE) {

		   genenamecol="Symbol", geneidcol=NULL, collaps=FALSE, 
		   llcol="CI.L", rlcol="CI.R", vcol=NULL, cvar=TRUE, 
		   metathr=0.01, jobname="MetaVolcano", outputfolder=".", 
		   draw='HTML', ncores=1, render = F) {

    if(!draw %in% c('PDF', 'HTML')) {

        stop("Oops! Seems like you did not provide a right 'draw' parameter.
              Try 'PDF' or 'HTML'")

    }

    # ---- Calculating variance from coifidence interval
    if(cvar == TRUE) {
        diffexp <- lapply(diffexp, function(...) calc_vi(..., llcol=llcol, rlcol=rlcol))
        vcol <- 'vi'
    } else {
        if(!is.null(llcol) && !is.null(rlcol)) {
            diffexp <- lapply(diffexp, function(...) calc_vi(..., llcol=llcol, rlcol=rlcol))
            vcol <- 'vi'
        } else {
            diffexp <- lapply(diffexp, function(...) calc_vi(..., llcol=NULL, rlcol=NULL,
                                                             vcol=vcol, foldchangecol=foldchangecol))
            vcol <- 'vi'
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
        dplyr::matches(paste(c(geneidcol, foldchangecol,
                               llcol, rlcol, vcol)[!sapply(c(geneidcol,
                                                             foldchangecol,
llcol, rlcol, vcol), is.null)], collapse = "|"))))

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

    remres <- do.call(rbind,
        mclapply(split(meta_diffexp, meta_diffexp[[genecol]]),
	    function(g) {
		remodel(g, foldchangecol, vcol)
                        }, mc.cores = ncores)
    )

    remres[[genecol]] <- rownames(remres)
    meta_diffexp <- merge(meta_diffexp, remres, by = genecol, all = TRUE)

    # --- Subsettig genes where REML doesnt converge
    meta_diffexp_err  <- dplyr::filter(meta_diffexp, error == TRUE)

    # --- Topconfects ranking
    meta_diffexp <- meta_diffexp %>%
        dplyr::filter(error != TRUE) # removing genes which REML
                                     # failed to converge

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

    # --- Keep all genes for the results report
    if(nrow(meta_diffexp_err) != 0) {
        meta_diffexp_err  <- meta_diffexp_err %>%
            dplyr::mutate(se=NA,
	                  index=NA,
		          `rank`=seq(nrow(meta_diffexp_err))+nrow(meta_diffexp))

        meta_diffexp  <- rbind(meta_diffexp, meta_diffexp_err)
    }

    meta_diffexp <- dplyr::arrange(meta_diffexp, `rank`)

    #####
    print(head(meta_diffexp))
    #####

    # --- Draw REM MetaVolcano
    gg <- plot_rem(meta_diffexp, jobname, outputfolder, genecol, metathr,
                      colors = colors, point_size = point_size,
                      label_genes = label_genes, label_top_n = label_top_n,
                      label_size = label_size, plot_title = plot_title,
                      show_legend = show_legend)

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

    if(render) {
      
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
		  featurefreq=ggplot()
		  )
    return(result)
}

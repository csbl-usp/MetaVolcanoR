#' A function to draw a forest plot from the REM MetaVolcano result
#'
#' This function draws a forest plot for a given gene based on the REM MetaVolcano result
#' @param gene query gene to plot
#' @param genecol name of the variable with genes <string>
#' @param remres data.table/data.frame output of the do.metafor function <data.table/data.frame>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param llcol left limit of the fold change coinfidence interval variable name <string>
#' @param rlcol right limit of the fold change coinfidence interval variable name <string>
#' @param studynames names of the input studies <character vector>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/ <string>
#' @param draw either 'PDF' or 'HTML' to save metaolcano as .pdf or .html respectively <string>
#' @keywords draw forest plot for a given gene
#' @export
#' @examples
#' draw.forest()
draw.forest <- function(gene, genecol, remres, foldchangecol, llcol, rlcol, studynames, jobname, outputfolder, draw) {
	remres %>%
		filter(!!rlang::sym(genecol) == gene) -> sremres
	
	stds <- unique(unlist(regmatches(colnames(sremres),
			regexec('_\\d+$', colnames(sremres)))))

	stds <- setNames(stds, studynames)
	
	edat <- Reduce(rbind, lapply(names(stds), function(sn) {
			std <- select(sremres, matches(paste0(genecol, '|', stds[sn])))
			colnames(std) <- gsub('_\\d+', '', colnames(std))
			std[['group']] <- sn
			std
		}))
	
	edat <- select(edat, c(!!rlang::sym(genecol), 
		       !!rlang::sym(foldchangecol), 
		       !!rlang::sym(llcol), 
		       !!rlang::sym(rlcol),
		       group))

	sdat <- data.frame(genecol = unique(edat[[genecol]]),
			   foldchangecol = sremres[['randomSummary']],
			   llcol = sremres[['randomCi.lb']],
			   rlcol = sremres[['randomCi.ub']],
			   group = 'FoldChange summary')
	
	colnames(sdat) <- c(genecol, foldchangecol, llcol, rlcol, 'group')
	dat <- rbind(edat, sdat)
	dat[['class']] <- ifelse(grepl('summary', dat[['group']]), "FoldChange summary", "Study")

	gg <- ggplot(dat, aes(x = group, y = !!rlang::sym(foldchangecol), 
		   color = `class`)) +
		geom_point() +
		geom_errorbar(aes(ymin = !!rlang::sym(llcol), 
				  ymax = !!rlang::sym(rlcol), 
			          width = 0.3,
				  color = `class`), alpha = 0.6) +
		scale_color_manual(values = c("#e6550d","#bdbdbd")) +
		scale_x_discrete(limits = rev(dat[['group']])) +
		theme_classic() +
		ggtitle(unique(edat[[genecol]])) +
		geom_hline(yintercept = 0, alpha=0.2, linetype = "dashed", size = 0.2) +
		theme(legend.position = "none") + 
		coord_flip()
	
	if(draw == "PDF") {

	  	pdf(paste0(normalizePath(outputfolder), 
		     "/Forestplot_", unique(edat[[genecol]]), '_', jobname, ".pdf"), width = 4, height = 5)
		plot(gg)
	  	dev.off()
  
	} else if (draw == "HTML") {

		htmlwidgets::saveWidget(as_widget(ggplotly(gg)), paste0(normalizePath(outputfolder),
							"/Forestplot_", unique(edat[[genecol]]), jobname, ".html"))
  	} else {
		stop("Seems like the draw parameter is invalid, try draw='PDF' or draw='HTML'")
	}	
}


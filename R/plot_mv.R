#' A MetaVolcano ploting function
#'
#' This function plots either the combining- or the vote-counting- MetaVolcanos
#' @param meta_diffexp data.frame/data.table containing the differential 
#'        expression inputs 
#' @param nstud the number of differential expression inputs <integer>
#' @param genecol column name of the variable to label genes in the .html 
#'        file <string>
#' @param comb wheather or not the drawing is for the combining-metavolcano 
#'        <logical>
#' @param metafc method for summarizing gene fold-changes across studies 
#'        c("Mean", "Median") <string>
#' @keywords draw metavolcano
#' @return \code{ggplot2} object
#' @export
#' @examples
#' plot_mv()
plot_mv <- function(meta_diffexp, nstud, genecol, comb, metafc) {
	if(comb) {
		
		# Drawing combining MetaVolcano
		g <- ggplot(meta_diffexp, aes(x = metafc, y = -log10(metap), 
					      text = !!rlang::sym(genecol))) +
        		geom_point(aes(color = degcomb), 
				   size = 0.5) +
        		labs(x = paste(metafc, "Fold Change"),
        		     y = "-log10(Fisher's(p values))")

        	} else {
	
		# Drawing vote-counting MetaVolcano
		g <- ggplot(meta_diffexp, aes(x = ddeg, y = ndeg, 
					      text = !!rlang::sym(genecol))) +
			geom_jitter(aes(color = degvcount), cex = 0.5, 
				    width = 0.45, height = 0.45) +
			scale_x_discrete(limits = -nstud:nstud) +
			labs(x = "Sign consistency",
			     y = "Number of times as differentially expressed")
	}

	g + theme_classic() +
		theme(panel.border= element_blank()) +
		theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) +
		theme(axis.line.x = element_line(color = "black", size = 0.6, 
						 lineend = "square"),
		      axis.line.y = element_line(color = "black", size = 0.6, 
						 lineend = "square")) +
		theme(legend.position = "none") +
		scale_color_manual(values=c("#377EB8", "grey", "#E41A1C"))
}

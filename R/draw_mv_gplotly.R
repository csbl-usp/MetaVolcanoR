#' A MetaVolcano drawing function
#'
#' This function creates a ggplotly object with the MetaVolcano
#' @param meta_geo2r data.frame/data.table containing all the inputed GEO2R outputs
#' @param nstud the number of inputed GEO2R outputs  <integer>
#' @param metathr the proportion of studies where a gene has to be consistently DE to be highlighted 'vote-counting meta-analysis' or the proportion of top meta-DEGs to be highlighted 'combining meta-analysis'<double>
#' @param genecol column name of the variable to label genes in the .html file <string>
#' @keywords draw metavolcano
#' @export
#' @examples
#' draw.mv.gplotly()
draw.mv.gplotly <- function(meta_geo2r, nstud, metathr, genecol) {
	ggplot(meta_geo2r, aes(x = metafc, y = -log10(metap), text = !!rlang::sym(genecol))) +
        	geom_point(aes(color = dcol_combin) , alpha = 0.9) +
        	labs(x = "Meta-Fold Change",
        	     y = "-log10(Meta-p value)") +
        	geom_hline(yintercept = quantile(as.numeric(-log10(meta_geo2r$metap)),
        	                                 metathr, na.rm = TRUE),
        	           linetype = "longdash", colour = "grey", alpha = 0.4) +
        	geom_vline(xintercept = c(-(quantile(as.numeric(meta_geo2r$metafc),
        	                                     1-metathr, na.rm = TRUE)),
        	                          quantile(as.numeric(meta_geo2r$metafc),
        	                                   1-metathr, na.rm = TRUE)),
        	           linetype = "longdash", colour = "grey", alpha = 0.4) +
		theme_classic() +
		theme(panel.border= element_blank()) +
		theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
		theme(axis.line.x = element_line(color = "black", size = 0.6, lineend = "square"),
		      axis.line.y = element_line(color = "black", size = 0.6, lineend = "square")) +
		theme(legend.position = "none") +
		scale_color_manual(values=c("#377EB8", "grey", "#E41A1C"))
}

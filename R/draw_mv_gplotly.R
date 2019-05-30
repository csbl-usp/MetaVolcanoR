#' A MetaVolcano drawing function
#'
#' This function creates a ggplotly object with the MetaVolcano
#' @param meta_geo2r data.frame/data.table containing all the inputed GEO2R outputs
#' @param nstud the number of inputed GEO2R outputs  <integer>
#' @param metathr top percentage of perturbed genes to be highlighted <double>
#' @param genecol column name of the variable to label genes in the .html file <string>
#' @param metap wheather or not the drawing is for the combining-metavolcano  <logical>
#' @param metafc method for summarizing gene fold-changes across studies c("Mean", "Median") <string>
#' @keywords draw metavolcano
#' @export
#' @examples
#' draw.mv.gplotly()
draw.mv.gplotly <- function(meta_geo2r, nstud, metathr, genecol, metap, metafc) {
	if(metap) {
		
		# Highlighting the top perturbed genes
		meta_geo2r <- meta_geo2r %>%
			mutate(idx = metafc*-log10(metap))
		meta_geo2r <- meta_geo2r %>%
			mutate(de_comb = ifelse(idx < quantile(meta_geo2r[['idx']], metathr/2), '0.Down-regulated',
						ifelse(idx > quantile(meta_geo2r[['idx']], 1-(metathr/2)), '2.Up-regulated', '1.Unperturbed')))
		# Drawing combining-MetaVolcano
		g <- ggplot(meta_geo2r, aes(x = metafc, y = -log10(metap), text = !!rlang::sym(genecol))) +
        		geom_point(aes(color = de_comb) , alpha = 0.9) +
        		labs(x = paste(metafc, "Fold Change"),
        		     y = "-log10(Fisher's(p values))")

        	} else {
		
		# Highlighting the top perturbed genes
		meta_geo2r <- meta_geo2r %>%
			mutate(idx = ddeg*ndeg)
		meta_geo2r <- meta_geo2r %>%
			mutate(de_comb = ifelse(idx < quantile(meta_geo2r[['idx']], metathr/2), '0.Down-regulated',
						ifelse(idx > quantile(meta_geo2r[['idx']], 1-(metathr/2)), '2.Up-regulated', '1.Unperturbed')))

		g <- ggplot(meta_geo2r, aes(x = ddeg, y = ndeg, text = !!rlang::sym(genecol))) +
			geom_jitter(aes(color = de_comb) , alpha = 0.9, width = 0.4) +
			labs(x = "Sign consistency",
			     y = "Number of times as differentially expressed")
	}

	g + theme_classic() +
		theme(panel.border= element_blank()) +
		theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
		theme(axis.line.x = element_line(color = "black", size = 0.6, lineend = "square"),
		      axis.line.y = element_line(color = "black", size = 0.6, lineend = "square")) +
		theme(legend.position = "none") +
		scale_color_manual(values=c("#377EB8", "grey", "#E41A1C"))
}

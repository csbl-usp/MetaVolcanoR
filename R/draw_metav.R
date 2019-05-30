#' A function to draw the 'Random Effect Model' MetaVolcano
#'
#' This function draws the 'Random Effect Model' MetaVolcano
#' @param meta_res data.frame/data.table containing all the inputed GEO2R outputs, and after REM calculation
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/ <string>
#' @param genecol column name of the variable to label genes in the .html file <string>
#' @param metathr top percentage of perturbed genes to be highlighted <double>
#' @param draw either 'PDF' or 'HTML' to save metaolcano as .pdf or .html respectively <string>
#' @keywords write 'combining meta-analysis' metavolcano
#' @export
#' @examples
#' draw.metav()
draw.metav <- function(meta_res, jobname, outputfolder, genecol, metathr, draw) {
  
  meta_res %>%
	  mutate(dircon2 = ifelse(`rank` <= quantile(meta_res[['rank']], metathr), dircon, NA)) %>%
	  filter(rank <  quantile(meta_res[['rank']], 0.6)) -> meta_res

  gg <- ggplot(arrange(meta_res, abs(randomSummary)),
	       aes(x = randomSummary, y = -log10(randomP), color = dircon2, text = !!rlang::sym(genecol))) +
	      geom_point() +
	      scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red", na.value = "grey80") +
	      labs(x = "Summary Fold Change",
	      	   y = "-log10(Summary p-value)",
		   color = "Sign consistency") +
	      geom_errorbarh(aes(xmax = randomCi.ub, xmin = randomCi.lb, color = dircon2), alpha = 0.6) +
	      theme_classic() +
	      theme(panel.border= element_blank()) +
	      theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
	      theme(axis.line.x = element_line(color = "black", size = 0.6, lineend = "square"),
	            axis.line.y = element_line(color = "black", size = 0.6, lineend = "square"))

  if(draw == "PDF") {

	  pdf(paste0(normalizePath(outputfolder), "/RandomEffectModel_MetaVolcano_", jobname, ".pdf"), width = 7, height = 10)
		  plot(gg)
	  dev.off()
  
  } else if (draw == "HTML") {

	  htmlwidgets::saveWidget(as_widget(ggplotly(gg)), paste0(normalizePath(outputfolder),
							"/RandomEffectModel_MetaVolcano_", jobname, ".html"))
  }	      
 }

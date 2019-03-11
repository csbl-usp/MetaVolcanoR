#' A function to draw the 'Random Effect Model' MetaVolcano 
#'
#' This function draws the 'Random Effect Model' MetaVolcano
#' @param meta_res data.frame/data.table containing all the inputed GEO2R outputs, and after REM calculation
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/ <string>
#' @keywords write 'combining meta-analysis' metavolcano
#' @export
#' @examples
#' draw.metavolcano.metap()
draw.metav <- function(meta_res, jobname, outputfolder) {
  gg <- ggplotly(
    ggplot(arrange(meta_res, abs(randomSummary)),
           aes(x = randomSummary, y = -log10(randomP), color = dircon, text = Gene.symbol)) +
      geom_point() +
      scale_color_gradient2(midpoint=0, low="blue", mid="white", high="red") +
      geom_errorbarh(aes(xmax = randomCi.ub, xmin = randomCi.lb, color = dircon)) +
      theme_classic() +
      theme(panel.border= element_blank()) +
      theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
      theme(axis.line.x = element_line(color = "black", size = 0.6, lineend = "square"),
            axis.line.y = element_line(color = "black", size = 0.6, lineend = "square"))
  )
  htmlwidgets::saveWidget(as_widget(gg), paste0(outputfolder, "RandomEffectModel_MetaVolcano_", jobname, ".html"))
}
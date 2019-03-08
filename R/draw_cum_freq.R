#' A function to visualize the inverse-cummulative DEG distribution
#'
#' This function create a ggplotly object with the inverse-cummulative DEG distribution
#' @param meta_geo2r data.frame/data.table containing all the inputed GEO2R outputs  
#' @param nstud the number of inputed GEO2R outputs  <integer>
#' @keywords draw inverse-cummulative DEG distribution
#' @export
#' @examples
#' draw.cum.freq()
draw.cum.freq <- function(meta_geo2r, nstud) {
  ggplotly(
    ggplot(cum.freq.data(meta_geo2r, nstud), aes(x = ndatasets, y = DEGs)) +
      geom_line(color = "lightblue", size = 2, alpha = 0.7) +
      theme_classic() +
      theme(panel.border= element_blank()) +
      theme(axis.text.x = element_text(angle=0, hjust = 1)) +
      theme(axis.line.x = element_line(color="black", size = 0.6, lineend = "square"),
            axis.line.y = element_line(color="black", size = 0.6, lineend = "square")) +
      guides(colour = guide_colorbar()) +
      labs(x = "Number of datasets",
           y = "Number of differentially expressed genes") +
      scale_x_discrete(limits=0:nstud)
  )
}
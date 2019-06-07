#' A function to visualize the inverse-cummulative DEG distribution
#'
#' This function create a ggplot object with the inverse-cummulative
#' DEG distribution
#' @param meta_diffexp data.frame/data.table containing all the input studies
#' @param nstud the number of inputed GEO2R outputs  <integer>
#' @keywords draw inverse-cummulative DEG distribution
#' @return \code{ggplot2} object
#' @export
#' @examples
#' draw_cum_freq()
draw_cum_freq <- function(meta_diffexp, nstud) {
    ggplot(cum_freq_data(meta_diffexp, nstud), aes(x = ndatasets, y = DEGs)) +
        geom_line(color = "#525252", size = 1) +
	geom_point(color = "#252525") +
        theme_classic() +
        theme(panel.border= element_blank()) +
        theme(axis.text.x = element_text(angle=0, vjust = 0.5)) +
        theme(axis.line.x = element_line(color="black", size = 0.6, 
					 lineend = "square"),
              axis.line.y = element_line(color="black", size = 0.6, 
					 lineend = "square")) +
        guides(colour = guide_colorbar()) +
        labs(x = "Number of datasets",
             y = "Number of differentially expressed genes") +
        scale_x_discrete(limits=0:nstud)
}

#' A function for DEG barplot visualization
#'
#' This function visualize as barplots the number of DEGs across the input 
#' studies
#' @param degbar_data output of the \code{set_degbar_data()} function 
#'        <data.fram/data.table>
#' @keywords draw DEG barplot
#' @return \code{ggplot2} object
#' @export
#' @examples
#' data(diffexplist)
#' diffexp <- lapply(diffexplist, function(...) deg_def(..., "pvalue", 
#'            "Log2FC", 0.05, 0))
#' bardat <- set_degbar_data(diffexp)
#' gg <- draw_degbar(bardat)
#' plot(gg)
draw_degbar <- function(degbar_data) {
    ggplot(degbar_data, aes(dataset)) +
        geom_bar(aes(fill = Regulation)) +
        theme_classic() +
        theme(panel.border= element_blank()) +
        theme(axis.text.x = element_text(angle=90, vjust = 0.5)) +
        theme(axis.line.x = element_line(color="black", size = 0.6, 
					 lineend = "square"),
              axis.line.y = element_line(color="black", size = 0.6, 
					 lineend = "square")) +
        guides(colour = guide_colorbar()) +
        labs(x = "Datasets",
             y = "Number of genes") +
        scale_fill_manual(values=c("#E41A1C", "grey", "#377EB8" )) +
        scale_x_discrete(labels=substr(unique(degbar_data[['dataset']]), 0, 30))
}

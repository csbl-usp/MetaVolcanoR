#' A function for feature barplot visualization
#'
#' This function visualize as barplots the number of features across the input 
#' studies
#' @param degbar_data output of the \code{set_degbar_data()} function 
#'        <data.fram/data.table>
#' @keywords draw feature barplot
#' @return \code{ggplot2} object
#' @export
#' @examples
#' data(diffexplist)
#' diffexp <- lapply(diffexplist, function(...) deg_def(..., "pvalue", 
#'            "Log2FC", 0.05, 0))
#' bardat <- set_degbar_data(diffexp)
#' gg <- draw_featurebar(bardat)
#' plot(gg)
draw_featurebar <- function(degbar_data, colors = c("#811820", "grey90", "#083e46")) {
    ggplot(degbar_data, aes(dataset)) +
        geom_bar(aes(fill = Regulation)) +
        theme_classic() +
        theme(panel.border= element_blank()) +
        theme(axis.text.x = element_text(angle=90, vjust = 0.5)) +
        theme(axis.line.x = element_line(color="black", linewidth = 0.6, 
					 lineend = "square"),
              axis.line.y = element_line(color="black", linewidth = 0.6, 
					 lineend = "square")) +
        guides(colour = guide_colorbar()) +
        labs(x = "Datasets",
             y = "Number of features") +
        scale_fill_manual(values=colors) +
        scale_x_discrete(labels=substr(unique(degbar_data[['dataset']]), 0, 30))
}

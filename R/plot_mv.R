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
#' @param colors vector of colors for the plot
#' @param point_size size of the points
#' @param label_genes character vector of specific genes to label
#' @param label_top_n number of top genes to label
#' @param label_size size of gene labels
#' @param plot_title custom plot title
#' @param show_legend whether to show the legend
#' @keywords draw metavolcano
#' @return \code{ggplot2} object
#' @export
#' @examples
#' data(diffexplist)
#' mv <- votecount_mv(diffexplist)
#' gg <- plot_mv(mv@metaresult, length(diffexplist), "Symbol", FALSE, "Mean")
#' plot(gg)
plot_mv <- function(meta_diffexp, nstud, genecol, comb, metafc,
                    colors = c("#083e46", "grey", "#811820"),
                    point_size = 0.5,
                    label_genes = NULL,
                    label_top_n = NULL,
                    label_size = 3,
                    plot_title = NULL,
                    show_legend = FALSE) {
    
    # Determine which genes to label
    genes_to_label <- c()
    if (!is.null(label_genes)) {
        genes_to_label <- union(genes_to_label, label_genes)
    }
    if (!is.null(label_top_n)) {
        if (comb) {
            top_genes <- meta_diffexp %>%
                dplyr::arrange(metap) %>%
                head(label_top_n) %>%
                dplyr::pull(!!rlang::sym(genecol))
        } else {
            top_genes <- meta_diffexp %>%
                dplyr::arrange(desc(abs(idx))) %>%
                head(label_top_n) %>%
                dplyr::pull(!!rlang::sym(genecol))
        }
        genes_to_label <- union(genes_to_label, top_genes)
    }
    
    # Add label column
    meta_diffexp <- meta_diffexp %>%
        dplyr::mutate(gene_label = ifelse(!!rlang::sym(genecol) %in% genes_to_label,
                                          as.character(!!rlang::sym(genecol)), 
                                          ""))
    
    if(comb) {
        # Drawing combining MetaVolcano
        g <- ggplot(meta_diffexp, aes(x = metafc, y = -log10(metap), 
                                      text = !!rlang::sym(genecol))) +
            geom_point(aes(color = degcomb), size = point_size) +
            labs(x = paste(metafc, "Fold Change"),
                 y = "-log10(Fisher's(p values))",
                 title = plot_title)
    } else {
        # Drawing vote-counting MetaVolcano
        g <- ggplot(meta_diffexp, aes(x = ddeg, y = ndeg, 
                                      text = !!rlang::sym(genecol))) +
            geom_jitter(aes(color = degvcount), size = point_size, 
                       width = 0.45, height = 0.45) +
            scale_x_continuous(breaks = -nstud:nstud, 
                   limits = c(-nstud - 0.5, nstud + 0.5)) +
            labs(x = "Sign consistency",
                 y = "Differential feature frequency",
                 title = plot_title)
    }
    
    g <- g + theme_classic() +
        theme(panel.border = element_blank()) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) +
        theme(axis.line.x = element_line(color = "black", linewidth = 0.6, 
                                        lineend = "square"),
              axis.line.y = element_line(color = "black", linewidth = 0.6, 
                                        lineend = "square")) +
        scale_color_manual(values = colors)
    
    # Add gene labels if specified
    if (length(genes_to_label) > 0) {
        g <- g + ggrepel::geom_text_repel(
            aes(label = gene_label),
            size = label_size,
            max.overlaps = Inf,
            box.padding = 0.5,
            point.padding = 0.3,
            segment.color = "grey50",
            show.legend = FALSE
        )
    }
    
    # Toggle legend
    if (!show_legend) {
        g <- g + theme(legend.position = "none")
    } else {
        g <- g + theme(legend.position = "right")
    }
    
    return(g)
}


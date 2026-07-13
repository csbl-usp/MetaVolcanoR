#' A function to plot the Random Effect Model (REM) MetaVolcano
#'
#' This function plots the REM MetaVolcano using ggplot2
#' @param meta_diffexp data.frame/data.table containing the REM results from 
#'        rem_mv() <data.table/data.frame>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/ <string>
#' @param genecol column name of the variable to label genes in the .html file
#'        <string>
#' @param metathr top percentage of perturbed genes to be highlighted <double>
#' @param colors named vector of colors for the gradient: c(low, mid, high, na)
#' @param point_size size of the points
#' @param label_genes character vector of specific genes to label
#' @param label_top_n number of top genes (by p-value) to label
#' @param label_size size of gene labels
#' @param plot_title custom plot title (NULL for no title)
#' @param show_legend whether to show the legend
#' @keywords write REM metavolcano
#' @return \code{ggplot2} object
#' @export
#' @examples
#' data(diffexplist)
#' diffexplist <- lapply(diffexplist, function(del) {
#'     dplyr::filter(del, grepl("MP", Symbol))
#' })
#' mv <- rem_mv(diffexplist, metathr = 0.1)
#' gg <- plot_rem(mv@metaresult, "MV", tempdir(), "Symbol", 0.01)
#' plot(gg)
plot_rem <- function(meta_diffexp, jobname, outputfolder, genecol, metathr,
                     colors = c(low = "#083e46", mid = "white", high = "#811820", na = "grey80"),
                     point_size = 0.6,
                     label_genes = NULL,
                     label_top_n = NULL,
                     label_size = 3,
                     plot_title = NULL,
                     show_legend = TRUE) {
    
    irank <- quantile(meta_diffexp[["rank"]], metathr)
    meta_diffexp %>%
        dplyr::mutate(signcon2 = ifelse(`rank` <= irank, signcon, NA)) %>%
        dplyr::mutate(Ci.ub = ifelse(`rank` <= irank, randomCi.ub, NA)) %>%
        dplyr::mutate(Ci.lb = ifelse(`rank` <= irank, randomCi.lb, NA)) %>%
        dplyr::filter(`rank` <  quantile(meta_diffexp[["rank"]], 0.6)) -> meta_res 
    
    # Determine which genes to label
    genes_to_label <- c()
    if (!is.null(label_genes)) {
        genes_to_label <- union(genes_to_label, label_genes)
    }
    if (!is.null(label_top_n)) {
        top_genes <- meta_res %>%
            dplyr::arrange(randomP) %>%
            head(label_top_n) %>%
            dplyr::pull(!!rlang::sym(genecol))
        genes_to_label <- union(genes_to_label, top_genes)
    }
    
    # Add label column
    meta_res <- meta_res %>%
        dplyr::mutate(gene_label = ifelse(!!rlang::sym(genecol) %in% genes_to_label,
                                          as.character(!!rlang::sym(genecol)), 
                                          ""))
    
    gg <- ggplot(dplyr::arrange(meta_res, abs(randomSummary)),
                 aes(x = randomSummary, y = -log10(randomP), color = signcon2, 
                     text = !!rlang::sym(genecol))) +
        geom_point(size = point_size) +
        scale_color_gradient2(midpoint = 0, 
                            low = colors["low"], 
                            mid = colors["mid"], 
                            high = colors["high"], 
                            na.value = colors["na"]) +
        labs(x = "Summary Fold-change",
             y = "-log10(Summary p-value)",
             color = "Sign consistency",
             title = plot_title) +
        geom_errorbarh(aes(xmax = Ci.ub, xmin = Ci.lb, color = signcon2)) +
        theme_classic() +
        theme(panel.border = element_blank()) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) +
        theme(axis.line.x = element_line(color = "black", linewidth = 0.6, 
                                        lineend = "square"),
              axis.line.y = element_line(color = "black", linewidth = 0.6, 
                                        lineend = "square"))
    
    # Add gene labels if specified
    if (length(genes_to_label) > 0) {
        gg <- gg + ggrepel::geom_text_repel(
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
        gg <- gg + theme(legend.position = "none")
    }
    
    return(gg)
}


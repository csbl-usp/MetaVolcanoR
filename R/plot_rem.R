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
#' @keywords write REM metavolcano
#' @export
#' @examples
#' plot_rem()
plot_rem <- function(meta_diffexp, jobname, outputfolder, genecol, metathr) {
    print(head(meta_diffexp))
    irank <- quantile(meta_diffexp[['rank']], metathr)
    meta_diffexp %>%
        dplyr::mutate(signcon2 = ifelse(`rank` <= irank, signcon, NA)) %>%
	dplyr::mutate(Ci.ub = ifelse(`rank` <= irank, randomCi.ub, NA)) %>%
	dplyr::mutate(Ci.lb = ifelse(`rank` <= irank, randomCi.lb, NA)) %>%
	dplyr::filter(`rank` <  quantile(meta_diffexp[['rank']], 0.6)) -> meta_res 


    ggplot(dplyr::arrange(meta_res, abs(randomSummary)),
        aes(x = randomSummary, y = -log10(randomP), color = signcon2, 
	    text = !!rlang::sym(genecol))) +
        geom_point() +
	scale_color_gradient2(midpoint=0, low="blue", mid="white", 
			      high="red", na.value = "grey80") +
	labs(x = "Summary Fold-change",
	     y = "-log10(Summary p-value)",
	     color = "Sign consistency") +
        geom_errorbarh(aes(xmax = Ci.ub, xmin = Ci.lb, 
			   color = signcon2), alpha = 0.6) +
        theme_classic() +
	theme(panel.border= element_blank()) +
	theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
	theme(axis.line.x = element_line(color = "black", size = 0.6, 
					 lineend = "square"),
	      axis.line.y = element_line(color = "black", size = 0.6, 
					 lineend = "square"))
}

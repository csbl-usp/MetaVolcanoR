#' @importFrom methods is
#' @importFrom stats setNames
NULL

#' A function to draw a forest plot from the REM MetaVolcano result
#'
#' This function draws a forest plot for a given gene based on the REM 
#' MetaVolcano result
#' @param remres MetaVolcano object. Output of the rem_mv() function
#'        <MetaVolcano>
#' @param gene query gene to plot
#' @param genecol name of the variable with genes <string>
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param llcol left limit of the fold change coinfidence interval variable
#'        name <string>
#' @param rlcol right limit of the fold change coinfidence interval variable
#'        name <string>
#' @param jobname name of the running job <string>
#' @param outputfolder /path where to write the results/ <string>
#' @param draw either 'PDF' or 'HTML' to save metaolcano as .pdf or .html
#'        respectively <string>
#' @param colors named vector of colors: c(positive, negative, neutral)
#' @param point_size size of the points
#' @param plot_width width of saved plot (inches for PDF, pixels for HTML)
#' @param plot_height height of saved plot (inches for PDF, pixels for HTML)
#' @param plot_title custom plot title (defaults to gene name)
#' @keywords draw forest-plot gene
#' @return \code{ggplot2} object
#' @export
#' @examples
#' data(diffexplist)
#' diffexplist <- lapply(diffexplist, function(del) {
#'     dplyr::filter(del, grepl("MP", Symbol))
#' })
#' mv <- rem_mv(diffexplist, metathr = 0.1)
#' gg <- draw_forest(mv, gene="MMP9")
#' plot(gg)
draw_forest <- function(remres, gene="MMP9", genecol="Symbol", 
                       foldchangecol="Log2FC", llcol="CI.L", rlcol="CI.R", 
                       jobname="MetaVolcano", outputfolder = tempdir(), draw="PDF",
                       colors = c(positive = "#811820", negative = "#083e46", 
                                 neutral = "#bdbdbd", reference = "#969696"),
                       point_size = 2,
                       plot_width = 4,
                       plot_height = 5,
                       plot_title = NULL) {
    
    if(!draw %in% c('PDF', 'HTML')) {
        stop("Oops! Seems like you did not provide a right 'draw' parameter. Try 'PDF' or 'HTML'")
    }
    
    if(is(remres) != "MetaVolcano") {
        stop("Oops! Please, provide a MetaVolcano object as input")
    }
    
    rem <- merge(remres@metaresult, remres@input, by = genecol) %>%
        dplyr::filter(!!rlang::sym(genecol) == gene) -> sremres
    
    if(nrow(sremres) == 0) {
        stop(paste("Oops! Seems that", gene, "is not in the provided REM result"))
    }
    
    stds <- unique(unlist(regmatches(colnames(sremres),
                                    regexec('_\\d+$', colnames(sremres)))))
    
    if(is.null(remres@inputnames)) {
        message("We recomend providing a character vector with the names of the input studies")
        stds <- setNames(stds, paste('study_', seq_along(stds)))
    } else {
        stds <- setNames(stds, remres@inputnames)
    }
    
    # setting data for visualization
    edat <- Reduce(rbind, lapply(names(stds), function(sn) {
        std <- dplyr::select(sremres, 
                           dplyr::matches(paste0(genecol, '|', stds[sn], '$')))
        colnames(std) <- gsub('_\\d+$', '', colnames(std))
        std[['group']] <- sn
        std
    }))
    
    if(!all(c(genecol, foldchangecol, llcol, rlcol) %in% colnames(edat))) {
        stop("Oops! Please, check the match among the provided parameters and the colnames of the remres@metaresult and remres@input")
    }
    
    edat <- dplyr::select(edat, c(!!rlang::sym(genecol), 
                                 !!rlang::sym(foldchangecol), 
                                 !!rlang::sym(llcol), 
                                 !!rlang::sym(rlcol),
                                 group))
    
    sdat <- data.frame(genecol = unique(edat[[genecol]]),
                      foldchangecol = sremres[['randomSummary']],
                      llcol = sremres[['randomCi.lb']],
                      rlcol = sremres[['randomCi.ub']],
                      group = 'FoldChange summary')
    
    colnames(sdat) <- c(genecol, foldchangecol, llcol, rlcol, 'group')
    dat <- rbind(edat, sdat)
    dat[['class']] <- ifelse(grepl('summary', dat[['group']]), 
                              "FoldChange summary", "Study")
    
    sumfc <- dplyr::filter(dat, grepl("summary", `class`))[[foldchangecol]]
    maxfc <- max(dat[[rlcol]])
    minfc <- min(dat[[llcol]])
    
    if(sumfc > 0) {
        sumcol <- unname(colors["positive"])
        minlim <- -maxfc
        maxlim <- maxfc
    } else {
        sumcol <- unname(colors["negative"])
        minlim <- minfc
        maxlim <- -minfc
    }
    
    # Use custom title or default to gene name
    plot_title_text <- if(is.null(plot_title)) unique(edat[[genecol]]) else plot_title
    
    gg <- ggplot(dat, aes(x = group, y = !!rlang::sym(foldchangecol), 
                         color = `class`)) +
        geom_point(size = point_size) +
        geom_errorbar(aes(ymin = !!rlang::sym(llcol), 
                         ymax = !!rlang::sym(rlcol), 
                         width = 0.1,
                         color = `class`)) +
        scale_color_manual(values = c("FoldChange summary" = sumcol, "Study" = unname(colors["neutral"]))) +
        scale_x_discrete(limits = rev(dat[['group']])) +
        theme_classic() +
        ggtitle(plot_title_text) +
        geom_hline(yintercept = 0, linetype = "solid", 
                  linewidth = 0.3, color = colors["reference"]) +
        geom_hline(yintercept = sumfc, linetype = "dashed", 
                  linewidth = 0.5, color = sumcol) +
        theme(legend.position = "none") + 
        scale_y_continuous(limits=c(minlim, maxlim)) +
        coord_flip()
    
    if(draw == "PDF") {  
        ggsave(filename = paste0(normalizePath(outputfolder), 
                                "/Forestplot_", unique(edat[[genecol]]), '_', jobname, ".pdf"),
               plot = gg,
               width = plot_width, 
               height = plot_height,
               device = "pdf")
    } else if (draw == "HTML") {
        htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
                               paste0(normalizePath(outputfolder),
                                     "/Forestplot_", unique(edat[[genecol]]),
                                     jobname, ".html"))
    }
    
    return(gg)
}


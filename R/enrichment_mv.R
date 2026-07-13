#' @importFrom stats setNames
NULL

#' Run gene set enrichment analysis on MetaVolcano results
#'
#' This function performs fast gene set enrichment analysis (fGSEA) on the
#' results of a random-effects model meta-analysis. Three ranking strategies
#' are supported: by summary fold-change, by signed significance, or by
#' a weighted combination of both. Gene sets can be provided directly or
#' downloaded automatically from MSigDB via msigdbr.
#' @param remres MetaVolcano object. Output of the rem_mv() function
#'        <MetaVolcano>
#' @param pathways named list of gene sets (character vectors). If NULL,
#'        gene sets are automatically downloaded via msigdbr using species,
#'        category, and subcategory parameters
#' @param species species name for msigdbr. Default "Homo sapiens"
#' @param category MSigDB category: "H" (hallmarks), "C2" (curated),
#'        "C5" (GO), "C7" (immunologic), etc. Default "C5"
#' @param subcategory optional MSigDB subcategory: "GO:BP", "GO:CC",
#'        "GO:MF", "CP:KEGG", "CP:REACTOME", etc. Default NULL (all)
#' @param ranking method for ranking genes: "fc" uses the REM summary
#'        fold-change; "signed_p" uses -log10(p) * sign(FC);
#'        "weighted_fc" uses FC * -log10(p) <string>
#' @param minSize minimum gene set size to test <integer>
#' @param maxSize maximum gene set size to test <integer>
#' @param nPermSimple number of permutations for p-value estimation <integer>
#' @param plot_top_n number of top pathways to plot <integer>
#' @param plot_padj adjusted p-value threshold for plotting. Default 0.1
#' @param clean_pathway_names logical. If TRUE, removes database prefixes
#'        (e.g. GOBP_, KEGG_) and replaces underscores with spaces.
#'        Default TRUE
#' @param colors vector of 2 colors for the plot c(down, up)
#' @param plot_title custom plot title. If NULL, auto-generated
#' @param seed random seed for reproducibility <integer>
#' @keywords enrichment fgsea gsea pathway
#' @return list with two elements: \code{result} (data.frame with fGSEA
#'         results) and \code{plot} (ggplot2 lollipop plot)
#' @export
#' @examples
#' \dontrun{
#' data(diffexplist)
#' mv <- rem_mv(diffexplist, metathr = 0.01)
#'
#' # Using automatic MSigDB download (GO Biological Process)
#' enrich <- enrichment_mv(mv, subcategory = "GO:BP")
#' enrich$plot
#'
#' # Using Hallmarks
#' enrich <- enrichment_mv(mv, category = "H")
#'
#' # Using KEGG
#' enrich <- enrichment_mv(mv, category = "C2", subcategory = "CP:KEGG")
#'
#' # Using custom pathways
#' library(msigdbr)
#' my_paths <- split(
#'     msigdbr(species = "Homo sapiens", category = "C5")$gene_symbol,
#'     msigdbr(species = "Homo sapiens", category = "C5")$gs_name
#' )
#' enrich <- enrichment_mv(mv, pathways = my_paths)
#' head(enrich$result)
#' }
enrichment_mv <- function(remres, 
                           pathways = NULL,
                           species = "Homo sapiens",
                           category = "C5",
                           subcategory = NULL,
                           ranking = "fc",
                           minSize = 15,
                           maxSize = 500,
                           nPermSimple = 10000,
                           plot_top_n = 30,
                           plot_padj = 0.1,
                           clean_pathway_names = TRUE,
                           colors = c(down = "#083e46", up = "#811820"),
                           plot_title = NULL,
                           seed = 42) {
    
    if (!requireNamespace("fgsea", quietly = TRUE)) {
        stop("Package 'fgsea' is required. Install with: BiocManager::install('fgsea')")
    }
    
    if (!methods::is(remres, "MetaVolcano")) {
        stop("Please provide a MetaVolcano object as input (output of rem_mv())")
    }
    
    meta <- remres@metaresult
    
    if (!all(c("randomSummary", "randomP", "Symbol") %in% colnames(meta))) {
        stop("Input must be a REM MetaVolcano result containing randomSummary, randomP, and Symbol columns")
    }
    
    # Build pathways from msigdbr if not provided
    if (is.null(pathways)) {
        if (!requireNamespace("msigdbr", quietly = TRUE)) {
            stop("Package 'msigdbr' is required when pathways = NULL. ",
                 "Install with: install.packages('msigdbr')")
        }
        
        msig_args <- list(species = species, category = category)
        if (!is.null(subcategory)) {
            msig_args$subcategory <- subcategory
        }
        
        msig_data <- do.call(msigdbr::msigdbr, msig_args)
        pathways <- split(msig_data$gene_symbol, msig_data$gs_name)
        
        message("Loaded ", length(pathways), " gene sets from MSigDB ",
                category, 
                ifelse(!is.null(subcategory), paste0(":", subcategory), ""))
    } else {
        if (!is.list(pathways) || is.null(names(pathways))) {
            stop("pathways must be a named list of character vectors")
        }
    }
    
    # Deduplicate before ranking (applies to all methods)
    meta_clean <- meta %>%
        dplyr::filter(!is.na(randomSummary), !is.na(randomP)) %>%
        dplyr::group_by(Symbol) %>%
        dplyr::slice_min(randomP, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()
    
    # Build gene ranking
    if (ranking == "fc") {
        gene_ranks <- setNames(meta_clean$randomSummary, meta_clean$Symbol)
    } else if (ranking == "signed_p") {
        gene_ranks <- setNames(
            -log10(meta_clean$randomP) * sign(meta_clean$randomSummary),
            meta_clean$Symbol
        )
    } else if (ranking == "weighted_fc") {
        gene_ranks <- setNames(
            meta_clean$randomSummary * -log10(meta_clean$randomP),
            meta_clean$Symbol
        )
    } else {
        stop("ranking must be 'fc', 'signed_p', or 'weighted_fc'")
    }
    
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)
    
    set.seed(seed)
    fgsea_res <- fgsea::fgsea(
        pathways    = pathways,
        stats       = gene_ranks,
        minSize     = minSize,
        maxSize     = maxSize,
        nPermSimple = nPermSimple
    ) %>%
        dplyr::arrange(pval) %>%
        dplyr::filter(!is.na(padj))
    
    if (is.null(plot_title)) {
        plot_title <- paste0("Functional enrichment (ranking: ", ranking, ")")
    }
    
    plot_df <- fgsea_res %>%
        dplyr::filter(padj < plot_padj) %>%
        dplyr::slice_min(padj, n = plot_top_n)
    
    if (nrow(plot_df) == 0) {
        message("No significant pathways found (padj < ", plot_padj, ")")
        return(list(result = fgsea_res, plot = ggplot2::ggplot()))
    }
    
    # Clean pathway names if requested
    if (clean_pathway_names) {
        plot_df <- plot_df %>%
            dplyr::mutate(
                pathway_clean = gsub("^[A-Z]+_", "", pathway),
                pathway_clean = gsub("_", " ", pathway_clean),
                pathway_clean = stringr::str_to_sentence(pathway_clean)
            )
    } else {
        plot_df <- plot_df %>%
            dplyr::mutate(pathway_clean = pathway)
    }
    
    plot_df <- plot_df %>%
        dplyr::mutate(
            direction     = ifelse(NES > 0, "Upregulated", "Downregulated"),
            pathway_clean = forcats::fct_reorder(pathway_clean, NES)
        )
    
    gg <- ggplot2::ggplot(plot_df,
                          ggplot2::aes(x = NES, y = pathway_clean, 
                                       colour = direction)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, 
                                            y = pathway_clean, 
                                            yend = pathway_clean),
                              linewidth = 0.7, alpha = 0.6) +
        ggplot2::geom_point(ggplot2::aes(size = -log10(padj)), alpha = 0.9) +
        ggplot2::scale_colour_manual(
            values = c("Upregulated"   = unname(colors["up"]),
                       "Downregulated" = unname(colors["down"])),
            name = "Direction"
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", 
                            colour = "grey40") +
        ggplot2::labs(x     = "NES",
                      y     = NULL,
                      title = plot_title,
                      size  = expression(-Log[10]~"(padj)")) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            panel.border = ggplot2::element_blank(),
            axis.line.x  = ggplot2::element_line(color = "black", 
                                                  linewidth = 0.6, 
                                                  lineend = "square"),
            axis.line.y  = ggplot2::element_line(color = "black", 
                                                  linewidth = 0.6, 
                                                  lineend = "square")
        )
    
    result_clean <- fgsea_res %>% dplyr::select(-leadingEdge)
    
    return(list(result = result_clean, plot = gg))
}
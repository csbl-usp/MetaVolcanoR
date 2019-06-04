#' A column renaming function merged inputs
#'
#' This function rename the columns of the merged inputs
#' @param diffexp list of data.frame/data.table (s) with DE results where lines
#'        are genes
#' @param genecol the column name of the geneID or gene name variable <string>
#' @param ncores the number of processors the user wants to use  <integer>
#' @keywords rename column
#' @export
#' @examples
#' rename_col()
rename_col <- function(diffexp, genecol, ncores) {
    ns <- names(diffexp)
    des <- mclapply(seq(diffexp), function(nstudy) {
        dex <- diffexp[[nstudy]]
        colnames(dex) <- paste(colnames(dex), nstudy, sep = "_")
        colnames(dex)[grep(genecol, colnames(dex))] <- genecol
        return(dex)
    }, mc.cores = ncores)
    names(des) <- ns
    return(des)
}

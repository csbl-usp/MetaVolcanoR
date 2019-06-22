#' A column renaming function merged inputs
#'
#' This function rename the columns of the merged inputs
#' @param diffexp list of data.frame/data.table (s) with DE results where lines
#'        are genes
#' @param genecol the column name of the geneID or gene name variable <string>
#' @keywords rename column
#' @return \code{data.tabledata.frame} with new colnames
#' @export
#' @examples
#' data(diffexplist)
#' lapply(diffexplist, colnames)
#' diffexp <- rename_col(diffexplist, "Symbol")
#' lapply(diffexp, colnames)
rename_col <- function(diffexp, genecol) {
    ns <- names(diffexp)
    des <- lapply(seq(diffexp), function(nstudy) {
        dex <- diffexp[[nstudy]]
        colnames(dex) <- paste(colnames(dex), nstudy, sep = "_")
        colnames(dex)[grep(genecol, colnames(dex))] <- genecol
        return(dex)
    })
    names(des) <- ns
    return(des)
}

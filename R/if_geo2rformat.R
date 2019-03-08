#' Checking if input files have GEO2R format
#'
#' Logical value
#' @param geo2r_res_file <GEO2R file>
#' @param inputfolder </input folder/>
#' @keywords input validation
#' @export
#' @examples
#' if.geo2rformat()
if.geo2rformat <- function(geo2r_res_file, inputfolder) {
  cnames <- gsub('"', '', unlist(strsplit(system(paste0("head -n 1 ", inputfolder, geo2r_res_file), intern = TRUE), '"\\t"')))
  if(length(which(cnames %in% c("ID", "adj.P.Val", "P.Value", "t", 
                                "B", "logFC", "Gene.symbol", "Gene.title", "CI.L", "CI.R"))) != 10) {
    message(paste(geo2r_res_file, "file does not look as a GEO2R output with Confidence Interval..."))
    FALSE
  } else {
    if(strsplit(system(paste0("wc ", inputfolder, geo2r_res_file), intern = TRUE), "\\s")[[1]][2] == 1) {
      message(paste(geo2r_res_file, "file looks as a GEO2R output but with no-rows"))
      FALSE
    } else {
      TRUE
    }
  }
}

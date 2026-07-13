# Register ggplot2 classes as S4-compatible
setOldClass(c("ggplot2::gg", "gg"))
setOldClass(c("ggplot2::ggplot", "ggplot", "ggplot2::gg", "gg"))

#' MetaVolcano S4 Class
#' 
#' @slot input data.frame
#' @slot inputnames character
#' @slot metaresult data.frame
#' @slot MetaVolcano gg
#' @slot featurefreq gg
#' @name MetaVolcano-class
#' @rdname MetaVolcano-class
#' @exportClass MetaVolcano
setClass("MetaVolcano", 
         slots = list(input = "data.frame",
                      inputnames = "character",
                      metaresult = "data.frame",
                      MetaVolcano = "gg",
                      featurefreq = "gg"))


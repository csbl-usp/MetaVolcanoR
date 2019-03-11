#' A function to model foldchange variance along several studies  
#'
#' This function calculate a Random Effect Model summary
#' @param gene named vector with foldchanges and variances (lines of the meta_geo2r: data.frame/data.table containing all the inputed GEO2R outputs)
#' @param foldchangecol the column name of the foldchange variable <string>
#' @keywords REM summary
#' @export
#' @examples
#' rnmodel()
rnmodel <- function(gene, foldchangecol) {
  fc <- gene[grep(foldchangecol, names(gene))]
  fc <- as.numeric(fc[which(!is.na(fc))])
  v <- gene[grep("vi", names(gene))]
  v <- as.numeric(v[which(!is.na(v))])
  
  # compute random and fixed effects from FCs e variances
  random <- tryCatch({ metafor::rma(fc, v, method="REML") }, 
                     error = function(e){ return(e) })
  # Increase iterations in case Fisher scoring algorithm doesn't converge
  if(any(class(random) == 'error')) {
    random <- tryCatch({ metafor::rma(fc, v, method="REML", 
                                      control = list(maxiter = 2000, stepadj = 0.5)) },
                       error = function(e){ print(e); return(e) })
  }
  # If metafor is still returning error, give up and register line for gene
  if(any(class(random) == 'error')) {
    df_res <- data.frame(dircon = length(which(fc > 0)) - length(which(fc < 0)),
                         randomSummary = NA, 
                         randomCi.lb = NA, 
                         randomCi.ub = NA, 
                         randomP = NA, 
                         het_QE = NA, 
                         het_QEp = NA,
                         het_QM = NA, 
                         het_QMp = NA,
                         #available = available, 
                         error = T)
  } else {
    df_res <- data.frame(dircon = length(which(fc > 0)) - length(which(fc < 0)),
                         randomSummary = as.numeric(random$beta), 
                         randomCi.lb = random$ci.lb, 
                         randomCi.ub = random$ci.ub, 
                         randomP = random$pval, 
                         het_QE = random$QE,
                         het_QEp = random$QEp,
                         het_QM = random$QM, 
                         het_QMp = random$QMp,
                         #available = available, 
                         error = F)
  }
  return(df_res)
}
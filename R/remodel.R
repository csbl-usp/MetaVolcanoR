#' A function to model foldchange variance along several studies  

#' This function calculate the REM-summary fold-change
#' @param gene named vector with foldchanges and variances <vector<
#' @param foldchangecol the column name of the foldchange variable <string>
#' @param vcol name of the fold change variance variable <string>
#' @keywords REM summary
#' @return \code{data.frame} with REM results for a gene
#' @export 
#' @examples
#' remodel()
remodel <- function(gene, foldchangecol, vcol) {

    fc <- as.numeric(dplyr::select(gene, matches(foldchangecol)))
    fc <- fc[which(!is.na(fc))]
    v <- as.numeric(dplyr::select(gene, matches(vcol)))
    v <- v[which(!is.na(v))]
    
    # compute random effect model for the fold-changes and its variace
    random <- tryCatch({ metafor::rma(fc, v, method="REML") }, 
                     error = function(e){ return(e) })

    # Increase iterations in case Fisher scoring algorithm doesn't converge
    if(any(is(random) == 'simpleError')) {
        random <- tryCatch({metafor::rma(fc, v, method="REML", 
                                      control = list(maxiter = 2000, 
						     stepadj = 0.5))},
                       error = function(e){ print(e); return(e) })
    }

    # If metafor is still returning error, give up and register line for gene
    if(any(is(random) == 'simpleError')) {
        df_res <- data.frame(signcon = length(which(fc>0))-length(which(fc<0)),
			 ntimes = length(fc != 0),
                         randomSummary = NA, 
                         randomCi.lb = NA, 
                         randomCi.ub = NA, 
                         randomP = NA, 
                         het_QE = NA, 
                         het_QEp = NA,
                         het_QM = NA, 
                         het_QMp = NA,
                         error = TRUE)
    } else {
        df_res <- data.frame(signcon = length(which(fc>0))-length(which(fc<0)),
			 ntimes = length(fc != 0),
			 randomSummary = as.numeric(random$beta), 
                         randomCi.lb = random$ci.lb, 
                         randomCi.ub = random$ci.ub, 
                         randomP = random$pval, 
                         het_QE = random$QE,
                         het_QEp = random$QEp,
                         het_QM = random$QM, 
                         het_QMp = random$QMp,
                         error = FALSE)
    }
    return(df_res)
}

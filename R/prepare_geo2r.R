#' Prepare a GEO2R result table for MetaVolcanoR
#'
#' Maps the columns exported by NCBI GEO2R (a limma-based differential
#' expression analysis run through the GEO web interface) to the format
#' expected by \code{rem_mv()}, \code{combining_mv()} and \code{votecount_mv()}.
#'
#' GEO2R tables do not include a confidence interval by default, but they do
#' report the moderated t-statistic. Because limma's moderated t is
#' \eqn{t = logFC / SE}, a standard error can be recovered as
#' \eqn{SE = logFC / t}, and a 95\% confidence interval as
#' \eqn{logFC \pm 1.96 \cdot SE}. This lets a plain GEO2R export feed the
#' random-effects model (REM) without any manual reformatting.
#'
#' @param diffexp A data.frame/data.table from GEO2R ("Download full table").
#' @param genenamecol Column holding the feature identifier (default
#'   "Gene.symbol"). GEO2R multi-symbol entries ("A///B") are split to the first.
#' @param logfccol Log2 fold-change column (default "logFC").
#' @param pvaluecol Raw p-value column (default "P.Value").
#' @param tcol Moderated t-statistic column used to derive the SE/CI
#'   (default "t"). Set to NULL to skip CI derivation (REM will then need a
#'   variance supplied another way).
#' @param collapse One of "min_pvalue" (keep, per symbol, the row with the
#'   smallest p-value; default) or "none" (leave duplicates untouched).
#' @param drop_unmapped Drop rows with missing/empty identifiers (default TRUE).
#'
#' @return A data.frame with columns: Symbol, Log2FC, pvalue, and (when tcol is
#'   available) CI.L, CI.R. Pass these to MetaVolcanoR via
#'   genenamecol="Symbol", foldchangecol="Log2FC", pcriteria="pvalue",
#'   llcol="CI.L", rlcol="CI.R".
#' @export
#'
#' @examples
#' \dontrun{
#'   tab <- read.delim("GSE12345.top.table.tsv", stringsAsFactors = FALSE)
#'   prep <- prepare_geo2r(tab)
#'   head(prep)
#' }
prepare_geo2r <- function(diffexp,
                          genenamecol  = "Gene.symbol",
                          logfccol     = "logFC",
                          pvaluecol    = "P.Value",
                          tcol         = "t",
                          collapse     = c("min_pvalue", "none"),
                          drop_unmapped = TRUE) {

  collapse <- match.arg(collapse)
  df <- as.data.frame(diffexp, stringsAsFactors = FALSE)

  needed <- c(genenamecol, logfccol, pvaluecol)
  missing <- setdiff(needed, colnames(df))
  if (length(missing)) {
    stop("prepare_geo2r(): columns not found in the GEO2R table: ",
         paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(df), collapse = ", "))
  }

  sym <- as.character(df[[genenamecol]])
  # GEO2R separates multiple mapped genes with "///" -> keep the first symbol.
  sym <- trimws(vapply(strsplit(sym, "///", fixed = TRUE),
                       function(x) if (length(x)) x[[1]] else NA_character_,
                       character(1)))

  out <- data.frame(
    Symbol = sym,
    Log2FC = suppressWarnings(as.numeric(df[[logfccol]])),
    pvalue = suppressWarnings(as.numeric(df[[pvaluecol]])),
    stringsAsFactors = FALSE
  )

  # Derive SE and 95% CI from the moderated t-statistic when available.
  if (!is.null(tcol) && tcol %in% colnames(df)) {
    tval <- suppressWarnings(as.numeric(df[[tcol]]))
    se <- ifelse(is.finite(tval) & tval != 0, out$Log2FC / tval, NA_real_)
    out$CI.L <- out$Log2FC - 1.96 * se
    out$CI.R <- out$Log2FC + 1.96 * se
  } else {
    message("prepare_geo2r(): no t-statistic column ('", tcol,
            "') found; CI.L/CI.R not derived. REM will require a variance ",
            "supplied by other means; Fisher and vote-counting are unaffected.")
  }

  if (drop_unmapped) {
    keep <- !is.na(out$Symbol) & nzchar(out$Symbol) &
            is.finite(out$Log2FC) & is.finite(out$pvalue)
    out <- out[keep, , drop = FALSE]
  }

  if (collapse == "min_pvalue" && anyDuplicated(out$Symbol)) {
    ord <- order(out$pvalue, na.last = NA)
    out <- out[ord, , drop = FALSE]
    out <- out[!duplicated(out$Symbol), , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

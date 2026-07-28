#' @importFrom utils head
#' @importFrom stats qnorm
NULL

# Declare global variables used in dplyr NSE
utils::globalVariables(c(
  "study", "value", "metap", "idx",
  "ndatasets", "DEGs", "dataset", "Regulation",
  "group", "degcomb", "ddeg", "ndeg", "degvcount",
  "signcon", "randomCi.ub", "randomCi.lb",
  "randomSummary", "randomP", "signcon2",
  "Ci.ub", "Ci.lb", "error", "index",
  # New variables from customization
  "gene_label",
  # Variables from prepare_* functions
  "log2FoldChange", "pvalue", "lfcSE", "Symbol", "Log2FC", "CI.L", "CI.R",
  "logFC", "PValue", "SE_approx",
  # Benjamini-Hochberg-adjusted p-values (rem_mv, combining_mv)
  "randomP.adjust", "metap.adjust",
  # Variables from enrichment_mv() (fgsea result columns, referenced via NSE)
  "pval", "padj", "pathway", "pathway_clean", "NES", "direction", "leadingEdge"
))


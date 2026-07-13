# MetaVolcanoR

**Meta-analysis visualization tool with publication-ready customization**

[![R](https://img.shields.io/badge/R-%3E%3D4.4.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-GPL--3-orange)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview

MetaVolcanoR combines differential gene, transcript, methylation and protein expression results from multiple studies to identify consistently perturbed genes. It implements three complementary meta-analysis strategies:

1. **Random Effects Model (REM)** - Rigorous statistical meta-analysis accounting for study variance
2. **Vote-counting** - Quick exploration of cross-study DEG consistency  
3. **Combining approach** - P-value aggregation using Fisher's method

All methods exploit volcano plot reasoning for intuitive visualization of meta-analysis results.

## Installation

### From GitHub (Latest Development Version)
```r
# Install devtools if needed
install.packages("devtools")

# Install MetaVolcanoR
devtools::install_github("iza-mcac/MetaVolcanoR")
```

### From Bioconductor (Coming Soon)
```r
# Will be available after Bioconductor submission
BiocManager::install("MetaVolcanoR")
```

## Quick Start
```r
library(MetaVolcanoR)

# Load example data (5 studies, ~20k genes each)
data(diffexplist)

# Run Random Effects Model meta-analysis
meta_results <- rem_mv(
  diffexp = diffexplist,
  metathr = 0.01,
  outputfolder = tempdir(),
  draw = "HTML"
)

# View interactive volcano plot
meta_results@MetaVolcano

# Explore forest plot for specific gene
draw_forest(meta_results, gene = "MMP9", draw = "PDF")
```

## New Features 🎨

**Plot Customization:**
- Custom color schemes (including colorblind-friendly palettes)
- Automatic or manual gene labeling
- Adjustable point sizes and plot dimensions
- Custom titles and legends

**Data Preparation Helpers:**
- `prepare_deseq2()` - One-line conversion from DESeq2 results
- `prepare_limma()` - Format limma/voom output
- `prepare_edger()` - Convert edgeR results
- `prepare_swish()` - Convert swish/Fishpond transcript-level results
- `prepare_limma()` - One-line conversion from limma/voom results
- `prepare_edger()` - One-line conversion from edgeR results

**Example with customization:**
```r
meta_custom <- rem_mv(
  diffexp = diffexplist,
  metathr = 0.01,
  # Customization options:
  colors = c(low = "navy", mid = "white", high = "darkred", na = "gray80"),
  label_genes = c("MMP9", "COL6A6", "MXRA5"),  # Specific genes
  label_size = 4,
  plot_title = "Disease vs Control Meta-Analysis",
  show_legend = TRUE
)
```

## Documentation

📖 **[Full Tutorial and Examples](https://iza-mcac.github.io/MetaVolcanoR/)** - Comprehensive vignette with:
- Data preparation from DESeq2/limma/edgeR
- All three meta-analysis methods explained
- Customization gallery
- Publication tips

## Input Data Requirements

Provide a **named list** of data frames, each containing:
- Gene identifiers (names or IDs)
- Log2 fold changes
- P-values
- *Optional*: Confidence intervals or variance (required for REM)

**Quick data prep from DESeq2:**
```r
library(DESeq2)
dds <- DESeq(dds)
res <- results(dds)

# One-line conversion!
deg_table <- prepare_deseq2(res)
```

## Three Meta-Analysis Methods

### 1. Random Effects Model (REM)

Most rigorous approach - accounts for between-study heterogeneity.
```r
rem_results <- rem_mv(diffexp = study_list, metathr = 0.01)
```

**When to use:** You have confidence intervals or standard errors

### 2. Vote-Counting

Fast exploration of cross-study DEG consistency.
```r
vote_results <- votecount_mv(
  diffexp = study_list,
  pvalue = 0.05,
  foldchange = 0.5
)
```

**When to use:** Quick overview, studies use different platforms

### 3. Combining Approach

Aggregates p-values (Fisher's method) and averages fold changes.
```r
comb_results <- combining_mv(
  diffexp = study_list,
  metafc = "Mean"  # or "Median"
)
```

**When to use:** Focus on statistical evidence aggregation


## License

GPL-3

---

**Contributors:** Izabela Mamede, Bárbara Marinho, Cesar Prada, Diogenes Lima, Helder Nakaya  
**Maintainer (at the moment):** Izabela Mamede (iza.mamede@gmail.com)

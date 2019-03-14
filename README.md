# MetaVolcanoR
Gene expression meta-analysis visualization tool

## Installation
### Clone repository
```
git clone https://github.com/csbl-usp/MetaVolcanoR
```
### In R
```
setwd('</parent folder where MetaVolcanoR was cloned/>')
devtools::install('MetaVolcanoR')
library(MetaVolcanoR)
```

## Running MetaVolcanoR


### Setting MetaVolcanoR parameters
```
inputfolder <- "/<path to input DE files>/"
outputfolder <- "/<path where outputs should be written>/"
jobname <- "my job"

pcriteria <- "pval column name"
genenamecol <- "gene name column name"
geneidcol <- "gene ID/probe/oligo/transcript column name"
foldchangecol <- "fold change column name"

pvalue <- 0.05 # pval threshold
logfc <- 0.0 # Fold change threshold

ncores <- 1 # number of processor user wants to use
collaps <- TRUE # c(TRUE, FALSE)

draw <- TRUE # weather or not to write the .html visualization
```

### Data input
```
infiles <- list.files(path = inputfolder)

infiles <- setNames(infiles, gsub("\\..+", "", infiles))

geo2r_res <- mclapply(infiles, function(f) {
  fread(paste0(inputfolder, f))
}, mc.cores = ncores)

```

### Draw DEGs by study and negative cumulative DEG distribution
```
ndegs <- draw.degbar.cum(geo2r_res, pcriteria, foldchangecol, genenamecol, geneidcol, pvalue, logfc, collaps, jobname, outputfolder, draw, ncores)
```

### Draw meta-volcano "Vote-counting approach"
```
metathr <- 0.8 # the proportion of studies/datasets/comparisons that a gene has to pass the pvalue and logfc thresholds to be highlighted <double>

meta_degs <- draw.metavolcano(geo2r_res, pcriteria, foldchangecol, genenamecol, geneidcol, pvalue, logfc, metathr, collaps, jobname, outputfolder, draw, ncores)
```

### Draw meta-volcano "Combining approach" fold-change Mean or Median & Fisher combining p-values
```
metap <- "Fisher"
metafc <- "Mean" # c("Mean", "Median")
metathr <- 0.8 # percentage of the top significant (smallest P-vals) and perturbed (extreme fold-changes) genes to be highlighted <double>

meta_degs_metap <- draw.metavolcano.metap(geo2r_res, pcriteria, foldchangecol, genenamecol, geneidcol, pvalue, logfc, metap, metafc, metathr, collaps, jobname, outputfolder, draw, ncores)

```

### Draw meta-volcano "Random effect model approach" calculating proper meta-Fold change
```
llcol <- "left limit of the fold change confidence interval variable name" # <string>
rlcol <- "right limit of the fold change confidence interval variable name" # <string>
vcol <- "name of the fold change variance variable" # <string>
cvar <- TRUE # weather or not to calculate the gene fold-change variance from the confidence interval limits <logical>

meta_degs_metafor <- do.metafor(geo2r_res, pcriteria, foldchangecol, genenamecol, geneidcol, pvalue, logfc, collaps, llcol, rlcol, vcol, cvar, jobname, outputfolder, draw, ncores)
```

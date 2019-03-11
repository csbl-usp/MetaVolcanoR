# MetaVolcanoR
Gene expression meta-analysis visualization tool

## Installation
### Clone repository
```
git clone csbl-usp/MetaVolcanoR
```
### In R
```
setwd('</parent folder where MetaVolcanoR was cloned/>')
devtools::install('MetaVolcanoR')
```

## Running MetaVolcanoR

### Required libraries
```
library(MetaVolcanoR)
library(data.table)
library(dplyr)
library(plotly)
```

### Setting MetaVolcanoR parameters
```
inputfolder <- "/<path to input DE files>/"
outputfolder <- "/<path where outputs should be written>/"
jobname <- "my job"
pcriteria <- "pval column name"
genenamecol <- "gene name column name"
geneidcol <- "gene ID/probe/oligo/transcript column name"
foldchangecol <- "fold change column name"
metathr <- 0.8 # percentage of DE to be considered as cDEG
pvalue <- 0.05 # pval threshold
logfc <- 0.0 # Fold change threshold
ncores <- 1 # number of processor user wants to use
collaps <- TRUE # c(TRUE, FALSE)
metap <- "Fisher"
metafc <- "Mean" # c("Mean", "Median")
cvar <- TRUE # c(TRUE, FALSE)
```

### data input
```
geo2r_res_files <- list.files(path = inputfolder)
geo2r_res_files <- setNames(geo2r_res_files, gsub("\\..+", "", geo2r_res_files))
geo2r_res <- mclapply(geo2r_res_files, function(...) fread(paste0(inputfolder, ...)),
                      mc.cores = ncores)
nstud <- length(geo2r_res)
```

### Draw DEGs by study and negative cumulative DEG distribution
```
meta_geo2r <- draw.degbar.cum(geo2r_res, pcriteria, foldchangecol, genenamecol, pvalue, logfc, collaps, jobname, outputfolder, ncores)
write.table(meta_geo2r, paste0(outputfolder, "deg_by_study_", jobname, ".tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
```

### Draw meta-volcano "Vote-counting approach"
```
meta_degs <- draw.metavolcano(meta_geo2r, genenamecol, metathr, nstud, jobname, collaps, outputfolder)
```

### Draw meta-volcano "Combining approach" FC Mean or Median & Fisher combining p-values
```
meta_degs_metap <- draw.metavolcano.metap(meta_geo2r, pcriteria, genenamecol, foldchangecol, metathr, nstud,
                                          jobname, collaps, metap, metafc, outputfolder)
write.table(meta_degs_metap, paste0(outputfolder, "meta_degs_", jobname, "_combining.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
```

### Draw meta-volcano "Random effect model approach" calculating proper meta-Fold change
```
meta_degs_metafor <- do.metafor(geo2r_res, pcriteria, genenamecol, geneidcol, foldchangecol, pvalue, logfc, collaps, jobname, outputfolder, ncores, cvar)
write.table(meta_degs_metafor, paste0(outputfolder, "metafor_degs_", jobname, "_summarizing.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
```

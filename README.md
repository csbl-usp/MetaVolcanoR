# MetaVolcanoR
Gene expression meta-analysis visualization tool

```
devtools::install("csbl/MetaVolcanoR")
library(MetaVolcanoR)

# Input parameters
inputfolder <- "/<whole path>/"
outputfolder <- "/<whlole path>"
jobname <- "myjob"
metathr <- 0.8 # percentage of DE to be considered as cDEG
pcriteria <-"P.Value" #c("adj.P.Val", "P.Value")
pvalue <- 0.05
logfc <- 0.1
ncores <- 4 # number of processor one wants to use
collaps <- TRUE # c(TRUE, FALSE)
metap <- "Fisher"
metafc <- "Mean" # c("Mean", "Median")
cvar <- TRUE # c(TRUE, FALSE)

# --- data input
geo2r_res_files <- list.files(path = inputfolder)
print(inputfolder)
geo2r_res_files <- setNames(geo2r_res_files, gsub("\\..+", "", geo2r_res_files))
message(head(geo2r_res_files))
geo2r_res_files <- geo2r_res_files[which(sapply(geo2r_res_files,
						function(...) if.geo2rformat(..., inputfolder)))] # checking files' format before read them into the R environment (Linux specific)
# reading input files      
geo2r_res <- mclapply(geo2r_res_files, function(...) fread(paste0(inputfolder, ...)),
		      mc.cores = ncores)

nstud <- length(geo2r_res)

# --- Running
# --- draw meta-volcano (stage 1) DEGs by study and cummulative DEG distribution

meta_geo2r <- draw.degbar.cum(geo2r_res, pcriteria, pvalue, logfc, collaps, jobname, outputfolder, ncores)
	# showing DEGs table
	print(head(meta_geo2r, 3))
	# writing table
	write.table(meta_geo2r, paste0(outputfolder, '/', "deg_by_study_", jobname, ".tsv"),
		    sep = "\t", row.names = FALSE, quote = FALSE)

# --- draw meta-volcano (stage 2) "Vote-counting aproach"
meta_degs <- draw.metavolcano(meta_geo2r, metathr, nstud, jobname, collaps, outputfolder)
	# showing meta-DEGs
	print(head(meta_degs, 3))

# --- draw meta-volcano (stage 2) "Combining aproach" FC Mean or Median & Fisher combining p-values
# --- [[[ here would be perfect to have a different slider for the metathr ]]]
meta_degs_metap <- draw.metavolcano.metap(meta_geo2r, pcriteria, metathr, nstud,
					  jobname, collaps, metap, metafc, outputfolder)
	# showing meta-DEGs
	print(head(meta_degs_metap, 3))
	# writing table
	write.table(meta_degs_metap, paste0(outputfolder, '/', "meta_degs_", jobname, "_combining.tsv"),
		    sep = "\t", row.names = FALSE, quote = FALSE)

# --- draw meta-volcano (stage 3) "Random effect model approach" calculating proper meta-Fold change
meta_degs_metafor <- do.metafor(geo2r_res, pcriteria, pvalue, logfc, collaps, jobname, outputfolder, ncores)
	print(head(meta_degs_metafor, 3))
	# writing table
	write.table(meta_degs_metafor, paste0(outputfolder, '/', "metafor_degs_", jobname, "_summarizing.tsv"),
		    sep = "\t", row.names = FALSE, quote = FALSE)
```

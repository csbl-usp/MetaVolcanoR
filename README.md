# MetaVolcanoR

Gene expression meta-analysis visualization tool.

## Overview

The MetaVolcanoR R package combines differential gene expression results. 
It implements three strategies to summarize gene expression activities from 
different studies. i) a naive vote-counting approach, ii) a combining-approach and 
iii) Random Effects Model (REM) approach. In all cases, MetaVolcano exploits the 
Volcano plot reasoning to visualize the gene expression meta-analysis results. 

## Installation
```
devtools::install_github('csbl-usp/MetaVolcanoR')
library(MetaVolcanoR)
```

## Usage
```
library(MetaVolcanoR,
	data.table,
	dplyr,
	parallel,
	plotly,
	metap,
	metafor,
	topconfects) 
```

### Input Data
Users should provide a named list of data.table/data.frame objects containing 
differential gene expression results. Each object of the list must contain the *gene name*,
*fold change*, and *p-value*. It is highly recomended to also include the *variance* or the
confidence intervals of the *fold change*. 

Take a look at the demo data. It includes differential gene expression results from five studies. 

```
data(diffexplist)
str(diffexplist)
head(diffexplist[[1]])
```

### Vote-counting approach
MetaVolcano identifies differential expressed genes (DEG) for each study based on the 
user-defined *p-value* and *fold change* thresholds. It displays the number of differentially
expressed and unperturbed genes per study. In addition, it plots the inverse cumulative distribution 
of the consistently DEG so the user can identify the number of genes whose expression is perturbed 
in at least 1 or n studies.

```
ndegs <- draw.degbar.cum(geo2r_res=diffexplist, 
			 pcriteria="pvalue",
			 foldchangecol="Log2FC",
			 genenamecol="Symbol", 
			 geneidcol="Symbol", 
			 pvalue=0.05, 
			 logfc=0, 
			 collaps=FALSE, 
			 jobname="demodata", 
			 outputfolder=".", 
			 draw="HTML", 
			 ncores=1)
```
![DEG by study and cummulative inverse distribution](https://.../.png)

MetaVolcano visualizes genes based on the number of studies where genes were identified as differentially
expressed and the their fold change *sign consistency*. It means that a gene that was differentially expressed 
in five studies and in three of them it was downregulated will get a *sign consistency* score of 2+(-3)=-1.
Based on the user selection, MetaVolcano can highligths the top *metathr* percentage of consistently perturbed genes.

```
meta_degs_vote <- draw.metavolcano(geo2r_res=diffexplist, 
				   pcriteria='pvalue',
				   foldchangecol='Log2FC',
				   genenamecol='Symbol',
				   geneidcol='Symbol',
				   pvalue=0.05,
				   logfc=0,
				   metathr=0.01,
				   collaps=FALSE, 
				   jobname="demodata",
				   outputfolder=".",
				   draw='HTML',
				   ncores=ncores)
```
![Vote-counting MetaVolcano](https://.../.png)

### Combining-approach 

The *combinig* MetaVolcano can also summarizes the *fold change* of a gene in different studies by the *mean* or *median* based 
on the user preference. In addition, the *combinig* MetaVolcano summarizes the gene differential expression *p-values* 
by mean of the Fisher method. The *combining* MetaVolcano can highligths the top *metathr* percentage of consistently perturbed genes.


```
meta_degs_comb <- draw.metavolcano.metap(geo2r_res=diffexplist, 
					 pcriteria='pvalue', 
					 foldchangecol='Log2FC',
					 genenamecol='Symbol', 
					 geneidcol='Symbol',
					 pvalue=0.05,
					 logfc=0, 
                                         metap='Fisher',
					 metafc='Mean',
					 metathr=0.01, 
					 collaps=FALSE,
					 jobname=jobname,
					 outputfolder=outputfolder,
					 draw='HTML',
					 ncores=1)
```
![Combining MetaVolcano](https://.../.png)


### Random Effect Model MetaVolcano

The *REM* MetaVolcano summarizes the gene fold change variance of several studies by a random effect model. 
Consequently, the summary p-value comes from the meta-analysis model and estimate the probability of the summary fold change
is equal zero. The *REM* MetaVolcano can highligths the top *metathr* percentage of consistently perturbed genes. This
perturbed ranking is done by the *topconfects* method. 


```
meta_degs_rem <- do.metafor(geo2r_res=diffexplist,
			    pcriteria='pvalue',
			    foldchangecol='Log2FC', 
			    genenamecol='Symbol',
			    geneidcol='Symbol',
			    pvalue=0.05,
			    logfc=0, 
			    collaps=FALSE,
			    llcol='CI.L',
			    rlcol='CI.R',
			    vcol=NULL, 
			    cvar=TRUE,
			    metathr=0.01,
			    jobname=jobname,
			    outputfolder=outputfolder, 
			    draw='HTML',
			    ncores=1)
```

![REM MetaVolcano](https://.../.png)
The *REM* MetaVolcano can also display the forest plot of a given gene based on the REM results.


```
draw.forest(gene="MXRA5", 
	    genecol="Symbol",
	    remres=meta_degs_rem,
	    foldchangecol="Log2FC",
	    llcol="CI.L",
	    rlcol="CI.R",
	    studynames=names(diffexplist),
	    jobname=jobname,
	    outputfolder=outputfolder,
	    draw="HTML")
```
![Forest plot](https://.../.png)


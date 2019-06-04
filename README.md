# MetaVolcanoR

Gene expression meta-analysis visualization tool.

## Overview

The MetaVolcanoR R package combines differential gene expression results. 
It implements three strategies to summarize gene expression activities from 
different studies. i) Random Effects Model (REM) approach. ii) a naive 
vote-counting approach, and iii) a combining-approach. In all cases, 
MetaVolcano exploits the Volcano plot reasoning to visualize the gene 
expression meta-analysis results.

## Installation
```
devtools::install_github('csbl-usp/MetaVolcanoR')
```

## Usage
Load required libraries.

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
differential gene expression results. Each object of the list must contain 
*gene name*, *fold change*, and *p-value* variables. It is highly recomended 
to also include *variance* or the confidence intervals of the *fold change* 
variables. 

Take a look at the demo data. It includes differential gene expression results
from five studies. 

```
data(diffexplist)
str(diffexplist)
head(diffexplist[[1]])
```

### Random Effect Model MetaVolcano

The *REM* MetaVolcano summarizes the gene fold change variance of several
studies by a random effect model. Consequently, the *summary p-value* comes
from the meta-analysis model and estimates the probability of the *summary 
fold-change* is equal zero. The *REM* MetaVolcano can highligths the top
*metathr* percentage of consistently perturbed genes. This perturbed ranking
is done by the *topconfects* (https://doi.org/10.1186/s13059-019-1674-7) method.


```
meta_degs_rem <- rem_mv(diffexp=diffexplist,
			pcriteria="pvalue",
			foldchangecol='Log2FC', 
			genenamecol='Symbol',
			geneidcol=NULL,
			collaps=FALSE,
			llcol='CI.L',
			rlcol='CI.R',
			vcol=NULL, 
			cvar=TRUE,
			metathr=0.01,
			jobname=jobname,
			outputfolder=outputfolder, 
			draw='HTML',
			ncores=ncores)
```
![REM MetaVolcano](https://github.com/csbl-usp/MetaVolcanoR/blob/master/REM_MV.png)

&nbsp;
The *REM* MetaVolcano can also display the forest plot of a given gene based 
on the REM results.

```
draw_forest(remres=meta_degs_rem,
	    gene="MMP9",
	    genecol="Symbol", 
	    foldchangecol="Log2FC",
	    llcol="CI.L", 
	    rlcol="CI.R",
	    studynames=names(diffexplist),
	    jobname=jobname,
	    outputfolder=outputfolder,
	    draw="HTML")

```
![Forest plot](https://github.com/csbl-usp/MetaVolcanoR/blob/master/forestplot.png)

### Vote-counting approach

MetaVolcano identifies differential expressed genes (DEG) for each study based 
on the user-defined *p-value* and *fold change* thresholds. It displays the 
number of differentially expressed and unperturbed genes per study. In addition,
it plots the inverse cumulative distribution of the consistently DEG so the
user can identify the number of genes whose expression is perturbed in at 
least 1 or n studies.

```
meta_degs_vote <- votecount_mv(diffexp=diffexplist,
			       pcriteria='pvalue',
			       foldchangecol='Log2FC',
			       genenamecol='Symbol',
			       geneidcol='Symbol',
			       pvalue=0.05,
			       foldchange=0, 
			       metathr=0.01,
			       collaps=FALSE,
			       jobname=jobname, 
			       outputfolder=outputfolder,
			       draw='HTML',
			       ncores=ncores)

```
![DEG by study and cummulative inverse distribution](https://github.com/csbl-usp/MetaVolcanoR/blob/dev/votecounting_pre_MV.png)

MetaVolcano visualizes genes based on the number of studies where genes were
identified as differentially expressed and the their fold change *sign 
consistency*. It means that a gene that was differentially expressed in five 
studies, from which three of them it was downregulated, will get a *sign 
consistency* score of *2 + (-3) = -1*. Based on the user selection, MetaVolcano
can highligths the top *metathr* percentage of consistently perturbed genes.

![Vote-counting MetaVolcano](https://github.com/csbl-usp/MetaVolcanoR/blob/master/votecounting_MV.png)

### Combining-approach 

The *combinig* MetaVolcano can also summarizes the *fold change* of a gene in
different studies by the *mean* or *median* based on the user preference. In
addition, the *combinig* MetaVolcano summarizes the gene differential
expression *p-values* by mean of the Fisher method. The *combining* MetaVolcano
can highligths the top *metathr* percentage of consistently perturbed genes.


```
meta_degs_comb <- combining_mv(diffexp=diffexplist,
			       pcriteria='pvalue', 
			       foldchangecol='Log2FC',
			       genenamecol='Symbol',
			       geneidcol='Symbol',
			       metafc='Mean',
			       metathr=0.01, 
			       collaps=TRUE,
			       jobname=jobname,
			       outputfolder=outputfolder,
			       draw='HTML',
			       ncores=ncores)

```
![Combining MetaVolcano](https://github.com/csbl-usp/MetaVolcanoR/blob/master/combining_MV.png)


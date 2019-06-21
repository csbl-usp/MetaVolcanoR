# MetaVolcanoR

Gene expression meta-analysis visualization tool.

## Overview

The MetaVolcanoR R package combines differential gene expression results. 
It implements three strategies to summarize gene expression activities from 
different studies. i) Random Effects Model (REM) approach. ii) a 
vote-counting approach, and iii) a combining-approach. MetaVolcano exploits 
the Volcano plot reasoning to visualize the gene expression meta-analysis 
results.

## Installation
```
BiocManager::install('MetaVolcanoR')
```

## Usage
Load required libraries.

```
library(MetaVolcanoR) 
```

### Input Data

Users should provide a named list of data.table/data.frame objects containing 
differential gene expression results. Each object of the list must contain 
*gene name*, *fold change*, and *p-value* variables. It is highly recomended 
to also include *variance* or the *confidence interval* of the *fold change* 
variables. 

Take a look at the demo data. It includes differential gene expression results
of five studies. 

```
data(diffexplist)
```

### Random Effect Model MetaVolcano

The *REM* MetaVolcano summarizes the gene fold change of several
studies taking into account the variance. The REM estimates a *summary p-value* 
which stand for the probability of the *summary fold-change* is not different
than zero. Users can set the *metathr* parameter to  highligth the top 
percentage of the most consistently perturbed genes. This perturbation 
ranking is defined following the  *topconfects* approach.


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
			jobname="MetaVolcano",
			outputfolder=".", 
			draw='HTML',
			ncores=5)

# REM results
head(meta_degs_rem@metaresult, 3)

# Ploting MetaVolcano
meta_degs_rem@MetaVolcano
```

&nbsp;
The *REM* MetaVolcano also allow users to explore the forest plot of a given 
gene based on the REM results.

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


```
draw_forest(remres=meta_degs_rem,
	    gene="COL6A6",
	    genecol="Symbol", 
	    foldchangecol="Log2FC",
	    llcol="CI.L", 
	    rlcol="CI.R",
	    studynames=names(diffexplist),
	    jobname=jobname,
	    outputfolder=outputfolder,
	    draw="HTML")

```


### Vote-counting approach

MetaVolcano identifies differential expressed genes (DEG) for each study based 
on the user-defined *p-value* and *fold change* thresholds. It displays the 
number of differentially expressed and unperturbed genes per study. In addition,
it plots the inverse cumulative distribution of the consistently DEG, so the
user can identify the number of genes whose expression is perturbed in at 
least 1 or n studies.

```
meta_degs_vote <- votecount_mv(diffexp=diffexplist,
			       pcriteria='pvalue',
			       foldchangecol='Log2FC',
			       genenamecol='Symbol',
			       geneidcol=NULL,
			       pvalue=0.05,
			       foldchange=0, 
			       metathr=0.01,
			       collaps=FALSE,
			       jobname="MetaVolcano", 
			       outputfolder=".",
			       draw='HTML',
			       ncores=4)

# Vote-counting results
head(meta_degs_vote@metaresult, 3)

# Plot DEG by study and DEG inverse cummulative distribution
meta_degs_vote@degfreq

```

The *vote-counting* MetaVolcano visualizes genes based on the number of studies 
where genes were identified as differentially expressed and the gene fold change
*sign consistency*. It means that a gene that was differentially expressed in 
five studies, from which three of them it was downregulated, will get a *sign 
consistency* score of *2 + (-3) = -1*. Based on the user preference, MetaVolcano
can highligths the top *metathr* percentage of consistently perturbed genes.

```
# Plot MetaVolcano
meta_degs_vote@MetaVolcano
```

### Combining-approach 

The *combinig* MetaVolcano summarizes the *fold change* of a gene in different
studies by the *mean* or *median* depending on the user preference. In addition, 
the *combinig* MetaVolcano summarizes the gene differential expression 
*p-values* using the Fisher method. The *combining* MetaVolcano can 
highligths the top *metathr* percentage of consistently perturbed genes.


```
meta_degs_comb <- combining_mv(diffexp=diffexplist,
			       pcriteria='pvalue', 
			       foldchangecol='Log2FC',
			       genenamecol='Symbol',
			       geneidcol=NULL,
			       metafc='Mean',
			       metathr=0.01, 
			       collaps=TRUE,
			       jobname="MetaVolcano",
			       outputfolder=".",
			       draw='HTML',
			       ncores=4)

# Combining results
head(meta_degs_comb@metaresult, 3)

# Plot MetaVolcano
meta_degs_comb@MetaVolcano

```


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapse_deg.R
\name{collapse_deg}
\alias{collapse_deg}
\title{A function to filter out geneIDs standing for the same gene name}
\usage{
collapse_deg(diffexp, genenamecol, pcriteria)
}
\arguments{
\item{diffexp}{data.frame/data.table output of the deg.def() function}

\item{genenamecol}{the column name of the gene name variable <string>}

\item{pcriteria}{the column name of the pvalue criteria to consider <string>}
}
\value{
\code{data.table} differential expression results with unique gene
        names
}
\description{
This function to remove redundant geneIDs standing for the same gene name
}
\examples{
data(diffexplist)
diffexp <- collapse_deg(diffexplist[[1]], "Symbol", "pvalue")
head(diffexp, 3)
}
\keyword{Probe}
\keyword{collapse}

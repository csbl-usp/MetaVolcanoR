% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_vi.R
\name{calc_vi}
\alias{calc_vi}
\title{A function to calculate variance from confidence interval limits}
\usage{
calc_vi(diffexp, llcol, rlcol)
}
\arguments{
\item{diffexp}{data.frame/data.table containing differential expression
results}

\item{llcol}{column name of the fold change coinfidence interval left limit
name <string>}

\item{rlcol}{column name of the fold change coinfidence interval left limit
name <string>}
}
\value{
\code{data.table/data.frame} with a new \code{vi} variable
}
\description{
This function takes the limits of a confidence interval (95%) 
a calculate a variance
}
\examples{
data(diffexplist)
diffexp <- calc_vi(diffexplist[[1]], "CI.L", "CI.R")
head(diffexp, 3)
}
\keyword{confidence}
\keyword{from}
\keyword{interval}
\keyword{variance}

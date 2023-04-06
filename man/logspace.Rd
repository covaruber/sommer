\name{logspace}
\alias{logspace}
\title{Decreasing logarithmic trend}
\description{
  logspace creates a vector with decreasing logaritmic trend.
}
\usage{
  logspace(n, start, end)
}

\arguments{
  \item{n}{number of values to generate.}
  \item{start}{initial value.}
  \item{end}{last value.}
}
\value{
  \describe{
    
    \item{$res}{ a vector of length n with logarithmic decrease trend.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

logspace(5, 1, .05)

}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmec}} }
\name{stackTrait}
\alias{stackTrait}
\title{Stacking traits in a dataset}
\description{
  \code{stackTrait} creates a dataset stacking traits in the long format to be used with the \code{\link{mmec}} solver for multi-trait models.
}
\usage{
  stackTrait(data, traits)
}

\arguments{
  \item{data}{a data frame with traits in wide format.}
  \item{traits}{variable names corresponding to the traits that should be in the long format.}
}
\value{
  \describe{
    
    \item{$res}{ a data frame with traits in long format.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_example)
DT <- DT_example
A <- A_example
head(DT)

DT2 <- stackTrait(DT, traits = c("Yield","Weight"))
head(DT2)

}
\seealso{The function \code{\link{vsc}} to know how to use \code{stackTrait} in the \code{\link{mmec}} solver.}
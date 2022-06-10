\name{isc}
\alias{isc}
\title{identity covariance structure}
\description{
  \code{isc} creates a diagonal covariance structure for the levels of the random effect  to be used with the \code{\link{mmec}} solver.
}
\usage{
  isc(x, thetaC=NULL, theta=NULL)
}

\arguments{
  \item{x}{vector of observations for the random effect.}
  \item{thetaC}{an optional matrix for constraints in the variance components.}
  \item{theta}{an optional matrix for initial values of the variance components.}
}
\value{
  \describe{
    \item{$res}{ a list with the provided vector and the variance covariance structure expected for the levels of the random effect.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
x <- as.factor(c(1:5,1:5,1:5));x
isc(x)
}
\seealso{See the function \code{\link{vsc}} to know how to use \code{isc} in the \code{\link{mmec}} solver.}
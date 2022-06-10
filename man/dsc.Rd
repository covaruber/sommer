\name{dsc}
\alias{dsc}
\title{diagonal covariance structure}
\description{
  \code{dsc} creates a diagonal covariance structure for the levels of the random effect  to be used with the \code{\link{mmec}} solver.
}
\usage{
  dsc(x, thetaC=NULL, theta=NULL)
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
dsc(x)
}
\seealso{See the function \code{\link{vsc}} to know how to use \code{dsc} in the \code{\link{mmec}} solver.}
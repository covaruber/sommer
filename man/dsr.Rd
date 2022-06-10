\name{dsr}
\alias{dsr}
\title{diagonal covariance structure}
\description{
  \code{dsr} creates a diagonal covariance structure for the levels of the random effect  to be used with the \code{\link{mmer}} solver.
}
\usage{
  dsr(x)
}

\arguments{
  \item{x}{vector of observations for the random effect.}
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
dsr(x)
}
\seealso{See the function \code{\link{vsr}} to know how to use \code{dsr} in the \code{\link{mmer}} solver.}
\name{usr}
\alias{usr}
\title{unstructured covariance structure}
\description{
  \code{usr} creates an unstructured covariance structure for specific levels of the random effect to be used with the \code{\link{mmer}} solver.
}
\usage{
  usr(x)
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
usr(x)
}
\seealso{The function \code{\link{vsr}} to know how to use \code{usr} in the \code{\link{mmer}} solver.}
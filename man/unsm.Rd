\name{unsm}
\alias{unsm}
\title{unstructured indication matrix}
\description{
  \code{unsm} creates a square matrix with ones in the diagonals and 2's in the off-diagonals to quickly specify an unstructured constraint in the Gtc argument of the \code{\link{vs}} function.
}
\usage{
  unsm(x)
}

\arguments{
  \item{x}{integer specifying the number of traits to be fitted for a given random effect.}
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
unsm(3)
}
\seealso{The function \code{\link{vs}} to know how to use \code{unsm} in the \code{\link{mmer}} solver.}
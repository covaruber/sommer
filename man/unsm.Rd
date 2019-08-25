\name{unsm}
\alias{unsm}
\title{unstructured indication matrix}
\description{
  \code{unsm} creates a square matrix with ones in the diagonals and 2's in the off-diagonals to quickly specify an unstructured constraint in the Gtc argument of the \code{\link{vs}} function.
}
\usage{
  unsm(x, reps=NULL)
}

\arguments{
  \item{x}{integer specifying the number of traits to be fitted for a given random effect.}
  \item{reps}{integer specifying the number of times the matrix should be repeated in a list format to provide easily the constraints in complex models that use the ds(), us() or cs() structures.}
}
\value{
  \describe{
    
    \item{$res}{ a matrix or a list of matrices with the constraints to be provided in the Gtc argument of the \code{\link{vs}} function.}
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
unsm(3,2)
}
\seealso{The function \code{\link{vs}} to know how to use \code{unsm} in the \code{\link{mmer}} solver.}
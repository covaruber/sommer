\name{uncm}
\alias{uncm}
\title{unconstrained indication matrix}
\description{
  \code{uncm} creates a square matrix with 2's in the diagnals and off-diagonals to quickly specify an unconstrained constraint in the Gtc argument of the \code{\link{vs}} function.
}
\usage{
  uncm(x, reps=NULL)
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
uncm(4)
}
\seealso{The function \code{\link{vs}} to know how to use \code{uncm} in the \code{\link{mmer}} solver.}
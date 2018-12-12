\name{list2usmat}
\alias{list2usmat}
\title{list or vector to unstructured matrix}
\description{
  \code{list2usmat} creates an unstructured square matrix taking a vector or list to fill the diagonal and upper triangular with the values provided.
}
\usage{
  list2usmat(sigmaL)
}

\arguments{
  \item{sigmaL}{vector or list of values to put on the matrix.}
}
\value{
  \describe{
    
    \item{$res}{ a matrix with the values provided.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
list2usmat(as.list(1:3))
list2usmat(as.list(1:10))
}
\seealso{The function \code{\link{vs}} to know how to use \code{list2usmat} in the \code{\link{mmer}} solver.}
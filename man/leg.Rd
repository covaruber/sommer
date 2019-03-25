\name{leg}
\alias{leg}
\title{Legendre polynomial matrix}
\description{
  Legendre polynomials of order 'n' are created given a vector 'x' and normalized to lay between values u and v.
}
\usage{
  leg(x,n=1,u=-1,v=1, intercept=TRUE, intercept1=FALSE)
}

\arguments{
  \item{x}{numeric vector to be used for the polynomial.}
  \item{n}{order of the Legendre polynomials.}
  \item{u}{lower bound for the polynomial.}
  \item{v}{upper bound for the polynomial.}
  \item{intercept}{a TRUE/FALSE value indicating if the intercept should be included.}
  \item{intercept1}{a TRUE/FALSE value indicating if the intercept should have value 1 (is multiplied by sqrt(2)).}
}
\value{
  \describe{
    
    \item{$S3}{ an Legendre polynomial matrix of order n.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
data("DT_legendre")
# you need to install the orthopolynom library
# leg(DT$X, n=1)
# leg(DT$X, n=2)

}
\seealso{The core functions of the package \code{\link{mmer}}}
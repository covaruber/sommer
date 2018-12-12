\name{ARMA}
\alias{ARMA}

\title{
Autocorrelation Moving average.
}
\description{
Creates an ARMA matrix of order one with parameters specified. 
}
\usage{
ARMA(x, rho=0.25, lambda=0.25)
}

\arguments{
  \item{x}{vector of the variable to define the factor levels for the ARMA covariance structure.}
  \item{rho}{rho value for the matrix.}
  \item{lambda}{dimensions of the square matrix.}
}
\details{
Specially useful for constructing covariance structures for rows and ranges to capture better the spatial variation trends in the field. The rho value is assumed fixed and values of the variance component will be optimized through REML.

}
\value{
If everything is defined correctly the function returns:
\describe{
\item{$nn}{the correlation matrix}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}

\examples{
x <- 1:4
R1 <- ARMA(x,rho=.25,lambda=0.2)
image(R1)
}
\seealso{The core function of the package \code{\link{mmer}}}

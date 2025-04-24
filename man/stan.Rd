\name{stan}
\alias{stan}

\title{
Standardize a vector of values in range 0 to 1
}
\description{
Simple function to map a vector of values to the range of 0 and 1 values to have a better behavior of the algorithm.
}
\usage{
stan(x, lb=0, ub=1)
}

\arguments{
  \item{x}{
A vector of numeric values.
}

  \item{lb}{
Lower bound value to map the x values.
}

  \item{ub}{
Upper bound value to map the x values.
}

}
\details{

Simple function to map a vector of values to the range of 0 and 1 values to have a better behavior of the algorithm.

}
\value{
\describe{
\item{$res}{new values in range 0 to 1}
}
}

\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

}

\examples{


x <- rnorm(20, 10, 3);x
stan(x)


}

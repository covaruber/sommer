\name{logspace}
\alias{logspace}
\title{Decreasing logarithmic trend}
\description{
  logspace creates a vector with decreasing logaritmic trend.
}
\usage{
  logspace(n, start=1, end=0.01, p =1, f=1)
}

\arguments{
  \item{n}{number of values to generate.}
  \item{start}{initial value.}
  \item{end}{last value.}
   \item{p}{power to be applied to the values.}
    \item{f}{a value from 1, 2, or 3 representing the different forms}
}
\value{
  \describe{
    
    \item{$res}{ a vector of length n with logarithmic decrease trend.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

plot(logspace(100,f=1))
plot(logspace(100,f=2))
plot(logspace(100,f=3))

}
\seealso{The core functions of the package \code{\link{mmes}}  }
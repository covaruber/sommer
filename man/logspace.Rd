\name{logspace}
\alias{logspace}
\title{Decreasing logarithmic trend}
\description{
  logspace creates a vector with decreasing logaritmic trend.
}
\usage{
  logspace(x, p=1)
}

\arguments{
  \item{x}{sequence of values to pass through the function.}
   \item{p}{power to be applied to the values.}
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

plot(logspace(1:100,p=1))
plot(logspace(1:100,p=2))
plot(logspace(1:100,p=3))

}
\seealso{The core functions of the package \code{\link{mmes}}  }
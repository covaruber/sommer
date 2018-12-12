\name{us}
\alias{us}
\title{unstructured covariance structure}
\description{
  \code{us} creates an unstructured covariance structure for specific levels of the random effect.
}
\usage{
  us(x)
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
data(DT_cpdata)
us(DT$Rowf)
}
\seealso{The function \code{\link{vs}} to know how to use \code{us} in the \code{\link{mmer}} solver.}
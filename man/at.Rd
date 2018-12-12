\name{at}
\alias{at}
\title{at covariance structure}
\description{
  \code{at} creates a diagonal covariance structure for specific levels of the random effect.
}
\usage{
  at(x, levs)
}

\arguments{
  \item{x}{vector of observations for the random effect.}
  \item{levs}{levels of the random effect to use for building the incidence matrices.}
}
\value{
  \describe{
    
    \item{$res}{ a list with the provided vector and the variance covariance structure expected.}
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
at(DT$Rowf)
at(DT$Rowf,c("1","2"))
}
\seealso{The function \code{\link{vs}} to know how to use \code{at} in the \code{\link{mmer}} solver. }
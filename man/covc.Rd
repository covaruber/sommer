\name{covc}
\alias{covc}
\title{covariance between random effects}
\description{
  \code{covc} merges the incidence matrices and covariance matrices of two random effects to fit an unstructured model between 2 different random effects to be fitted with the \code{\link{mmec}} solver.
}
\usage{
  covc(ran1, ran2)
}

\arguments{
  \item{ran1}{the random call of the first random effect.}
  \item{ran2}{the random call of the first random effect.}
}
\details{
This implementation aims to fit models where covariance between random variables is expected to exist. For example, indirect genetic effects.

}
\value{
  \describe{
    
    \item{$Z}{ a incidence matrix Z* = Z Gamma which is the original incidence matrix for the timevar multiplied by the loadings.}
  }
}
\references{

  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744
  
  Bijma, P. (2014). The quantitative genetics of indirect genetic effects: a selective review of modelling issues. Heredity, 112(1), 61-69.
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_ige)
DT <- DT_ige
covRes <- with(DT, covc( vsc(isc(focal)) , vsc(isc(neighbour)) ) )
str(covRes)
# look at DT_ige help page to see how to fit an actual model

}
\seealso{The function \code{\link{vsc}} to know how to use \code{covc} in the \code{\link{mmec}} solver.}
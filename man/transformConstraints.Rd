\name{transformConstraints}
\alias{transformConstraints}
\title{transformConstraints}
\description{
  \code{transformConstraints} takes a list of matrices with constraints and transforms all the non-zero values to the value desired. The purpose of this function is to make easy the transformation of initial constraints to a fixed-constraint list to be provided to a mixed model fitted with the mmer function.
}
\usage{
  transformConstraints(list0,value=1)
}

\arguments{
  \item{list0}{a list of matrices with constraints according to the rules specified in the vs function (0: not to be estimated, 1: positive, 2:unconstrained, 3:fixed).}
  \item{value}{value to be used to replace all the non-zero values in the constraint matrices.}
}
\value{
  \describe{
    
    \item{$res}{ a list with the modified constraint matrices.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
(a <- list(unsm(4), diag(4)))
transformConstraints(a, value=3)
}
\seealso{The function \code{\link{vs}} to know how to use \code{transformConstraints} in the \code{\link{mmer}} solver. }
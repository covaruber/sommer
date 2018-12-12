\name{pedtoK}
\alias{pedtoK}
\title{Pedigree to matrix}
\description{
  This function creates takes the inverse of the asreml Ainverse function (ginv element) and creates the additive relationship matrix to use it as a covariance matrix for a random effect. Other packages that allows you to obtain an additive relationship matrix from a pedigree is the `pedigreemm` package.
}
\usage{
  pedtoK(x, type="asreml")
}

\arguments{
  \item{x}{ginv element, output from the Ainverse function.}
  \item{type}{program used to obtain the pedigree. Currently only asreml.}
}
\value{
    
    \item{K}{ pedigree transformed in a relationship matrix.}
    \item{Kinv}{ inverse of the pedigree transformed in a relationship matrix.}
  
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{

# Ks <- pedtoK(asreml.model$ginv)
# A <- Ks$K

}
\seealso{The core functions of the package \code{\link{mmer}} }
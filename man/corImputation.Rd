\name{corImputation}
\alias{corImputation}
\title{Imputing a matrix using correlations}
\description{
  corImputation imputes missing data based on the correlation that exists between row levels.
}
\usage{
  corImputation(wide, Gu=NULL, nearest=1, roundR=FALSE)
}

\arguments{
  \item{wide}{numeric matrix with individuals in rows and time variable in columns (e.g., environments, genetic markers, etc.).}
  \item{Gu}{optional correlation matrix between the individuals or row levels. If NULL it will be computed as the correlation of t(wide).}
  \item{nearest}{integer value describing how many nearest neighbours (the ones showing the highest correlation) should be used to average and return the imputed value.}
  \item{roundR}{a TRUE/FALSE statement describing if the average result should be rounded or not. This may be specifically useful for categorical data in the form of numbers (e.g., -1,0,1).}
}
\value{
  \describe{
    
    \item{$res}{ a list with the imputed matrix and the original matrix.}
    
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{



}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmec}} }
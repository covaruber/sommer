\name{imputev}
\alias{imputev}
\title{Imputing a numeric or character vector}
\description{
This function is a very simple function to impute a numeric or character vector with the mean or median value of the vector.
}
\usage{
imputev(x, method="median")
}

\arguments{
  \item{x}{a numeric or character vector}
  \item{method}{the method to choose between mean or median}
}
\value{
\describe{

\item{$x}{a numeric or character vector imputed with the method selected.}

}
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\author{
Giovanny Covarrubias-Pazaran
}

\examples{
####=========================================####
#### generate your mickey mouse -log10(p-values)
####=========================================####
set.seed(1253)
x <- rnorm(100)
x[sample(1:100,10)] <- NA
imputev(x)
}
\seealso{The core function of the package \code{\link{mmer}}}
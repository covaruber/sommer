\name{dfToMatrix}
\alias{dfToMatrix}
\title{data frame to matrix}
\description{
  This function takes a matrix that is in data frame format and transforms it into a matrix. Other packages that allows you to obtain an additive relationship matrix from a pedigree is the `pedigreemm` package.
}
\usage{
  dfToMatrix(x, row="Row",column="Column",
             value="Ainverse", returnInverse=FALSE, 
             bend=1e-6)
}

\arguments{
  \item{x}{ginv element, output from the Ainverse function.}
  \item{row}{name of the column in x that indicates the row in the original relationship matrix.}
  \item{column}{name of the column in x that indicates the column in the original relationship matrix.}
  \item{value}{name of the column in x that indicates the value for a given row and column in the original relationship matrix.}
  \item{returnInverse}{a TRUE/FALSE value indicating if the inverse of the x matrix should be computed once the data frame x is converted into a matrix.}
    \item{bend}{a numeric value to add to the diagonal matrix in case matrix is singular for inversion.}
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

library(Matrix)
m <- matrix(1:9,3,3)
m <- tcrossprod(m)

mdf <- as.data.frame(as.table(m))
mdf

dfToMatrix(mdf, row = "Var1", column = "Var2", 
            value = "Freq",returnInverse=FALSE )


}
\seealso{The core functions of the package \code{\link{mmer}} }
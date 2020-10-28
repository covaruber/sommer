\name{dfToMatrix}
\alias{dfToMatrix}
\title{Pedigree to matrix}
\description{
  This function creates takes the inverse of the asreml Ainverse function (ginv element) and creates the additive relationship matrix to use it as a covariance matrix for a random effect. Other packages that allows you to obtain an additive relationship matrix from a pedigree is the `pedigreemm` package.
}
\usage{
  dfToMatrix(x, row="Row",column="Column",value="Ainverse", returnInverse=TRUE)
}

\arguments{
  \item{x}{ginv element, output from the Ainverse function.}
  \item{row}{name of the column in x that indicates the row in the original relationship matrix.}
  \item{column}{name of the column in x that indicates the column in the original relationship matrix.}
  \item{value}{name of the column in x that indicates the value for a given row and column in the original relationship matrix.}
  \item{returnInverse}{a TRUE/FALSE value indicating if the inverse of the x matrix should be computed once the data frame x is converted into a matrix.}
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

m <- matrix(1:9,3,3)
m <- m%*%t(m)

mdf <- as.data.frame(as.table(m))
mdf

dfToMatrix(mdf, row = "Var1", column = "Var2", 
            value = "Freq",returnInverse=F )


}
\seealso{The core functions of the package \code{\link{mmer}} }
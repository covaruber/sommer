\name{overlay}
\alias{overlay}
\title{Overlay Matrix}
\description{
  `overlay` adds r times the design matrix for model term t to the existing design matrix. Specifically, if the model up to this point has p effects and t has a effects, the a columns of the design matrix for t are multiplied by the scalar r (default value 1.0). This can be used to force a correlation of 1 between two terms as in a diallel analysis.
}
\usage{
  overlay(..., rlist=NULL, prefix=NULL, sparse=FALSE)
}

\arguments{
  \item{...}{as many vectors as desired to overlay.}
  \item{rlist}{a list of scalar values indicating the times that each incidence matrix overlayed should be multiplied by. By default r=1.}
  \item{prefix}{a character name to be added before the column names of the final overlay matrix. This may be useful if you have entries with names starting with numbers which programs such as asreml doesn't like, or for posterior extraction of parameters, that way 'grep'ing is easier.}
  \item{sparse}{a TRUE/FALSE statement specifying if the matrices should be built as sparse or regular matrices.}
}
\value{
  \describe{
    
    \item{$S3}{ an incidence matrix with as many columns levels in the vectors provided to build the incidence matrix.}
  }
}
\references{
  Fikret Isik. 2009. Analysis of Diallel Mating Designs. North Carolina State University, Raleigh, USA.
  
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
  
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
data("DT_halfdiallel")
DT <- DT_halfdiallel
head(DT)
DT$femalef <- as.factor(DT$female)
DT$malef <- as.factor(DT$male)
DT$genof <- as.factor(DT$geno)

with(DT, overlay(femalef,malef, sparse = TRUE))
with(DT, overlay(femalef,malef, sparse = FALSE))

}
\seealso{The core functions of the package \code{\link{mmer}} and a function for creating dummy variables for diallel models named \code{\link{add.diallel.vars}}. }
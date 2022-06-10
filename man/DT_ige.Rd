\name{DT_ige}
\alias{DT_ige}
\alias{A_ige}
\docType{data}
\title{
Data to fit indirect genetic effects.
}
\description{

This dataset contains phenotpic data for 98 individuals where they are measured with the purpose of identifying the effect of the neighbour in a focal individual. 

}
\usage{data("DT_ige")}
\format{
  The format is:
 chr "DT_ige"
}
\source{
This data was masked from a shared study.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
####=========================================####
#### EXAMPLES
#### Different models with sommer
####=========================================####

data(DT_ige)
DT <- DT_ige
Af <- A_ige
An <- A_ige
### Direct genetic effects model
# modDGE <- mmer(trait ~ block,
#                random = ~ focal,
#                rcov = ~ units, 
#                data = DT)
# summary(modDGE)$varcomp
# 
### Indirect genetic effects model without covariance between DGE and IGE
# modDGE <- mmer(trait ~ block,
#                random = ~focal + neighbour,
#                rcov = ~ units, 
#                data = DT)
# summary(modDGE)$varcomp
# 
### Indirect genetic effects model with covariance between DGE and IGE
# modIGE <- mmer(trait ~ block,
#                random = ~ gvsr(focal, neighbour),
#                rcov = ~ units, iters=4,
#               data = DT)
# summary(modIGE)$varcomp
#
### Indirect genetic effects model with covariance between DGE and IGE using relationship matrices
# modIGEb <- mmer(trait ~ block,
#                random = ~ gvsr(focal, neighbour, Gu=list(Af,An)),
#                rcov = ~ units, 
#                data = DT)
# summary(modIGEb)$varcomp
}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmec}} }

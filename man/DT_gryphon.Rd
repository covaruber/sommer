\name{DT_gryphon}
\alias{DT_gryphon}
\alias{A_gryphon}
\alias{P_gryphon}
\docType{data}
\title{
Gryphon data from the Journal of Animal Ecology
}
\description{
This is a dataset that was included in the Journal of animal ecology by Wilson et al. (2010; see references) to help users understand how to use mixed models with animal datasets with pedigree data. 

The dataset contains 3 elements:

gryphon; variables indicating the animal, the mother of the animal, sex of the animal, and two quantitative traits named 'BWT' and 'TARSUS'.

pedi; dataset with 2 columns indicating the sire and the dam of the animals contained in the gryphon dataset.

A; additive relationship matrix formed using the 'getA()' function used over the pedi dataframe.

}
\usage{data("DT_gryphon")}
\format{
  The format is:
 chr "DT_gryphon"
}
\source{
This data comes from the Journal of Animal Ecology. Please, if using this data cite Wilson et al. publication. If using our mixed model solver please cite Covarrubias' publication.
}
\references{

Wilson AJ, et al. (2010) An ecologist's guide to the animal model. Journal of Animal Ecology 79(1): 13-26.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
 
}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples using
#### command + shift + C |OR| control + shift + C
####=========================================####
# data(DT_gryphon)
# DT <- DT_gryphon
# A <- A_gryphon
# P <- P_gryphon
# #### look at the data
# head(DT)
# #### fit the model with no fixed effects (intercept only)
# mix1 <- mmer(BWT~1,
#              random=~vsr(ANIMAL,Gu=A),
#              rcov=~units,
#              data=DT)
# summary(mix1)$varcomp
# 
# ## mmec uses the inverse of the relationship matrix
# Ai <- as(solve(A + diag(1e-4,ncol(A),ncol(A))), Class="dgCMatrix")
# mix1b <- mmec(BWT~1,
#              random=~vsc(isc(ANIMAL),Gu=Ai),
#              rcov=~units, tolParConv = 1e-5,
#              data=DT)
# summary(mix1b)$varcomp
# 
# #### fit the multivariate model with no fixed effects (intercept only)
# mix2 <- mmer(cbind(BWT,TARSUS)~1,
#              random=~vsr(ANIMAL,Gu=A),
#              rcov=~vsr(units),
#              na.method.Y = "include2",
#              data=DT)
# summary(mix2)
# cov2cor(mix2$sigma$`u:ANIMAL`)
# cov2cor(mix2$sigma$`u:units`)

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmec}} }

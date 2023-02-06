\name{DT_cpdata}
\alias{DT_cpdata}
\alias{GT_cpdata}
\alias{MP_cpdata}
\docType{data}
\title{
Genotypic and Phenotypic data for a CP population
}
\description{
A CP population or F1 cross is the designation for a cross between 2 highly heterozygote individuals; i.e. humans, fruit crops, bredding populations in recurrent selection.

This dataset contains phenotpic data for 363 siblings for an F1 cross. These are averages over 2 environments evaluated for 4 traits; color, yield, fruit average weight, and firmness. The columns in the CPgeno file are the markers whereas the rows are the individuals. The CPpheno data frame contains the measurements for the 363 siblings, and as mentioned before are averages over 2 environments.

}
\usage{data("DT_cpdata")}
\format{
  The format is:
 chr "DT_cpdata"
}
\source{
This data was simulated for fruit breeding applications.
}
\references{

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

}
\examples{

# ####=========================================####
# #### For CRAN time limitations most lines in the
# #### examples are silenced with one '#' mark,
# #### remove them and run the examples using
# #### command + shift + C |OR| control + shift + C
# ####=========================================####
#
# data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix
# A <- A.mat(GT) # additive relationship matrix
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmer(Yield~1,
#               random=~vsr(id,Gu=A)
#                       + Rowf + Colf,
#               rcov=~units,
#               data=DT)
# summary(mix1)$varcomp
#
# ## mmec uses the inverse of the relationship matrix
# Ai <- as(solve(A + diag(1e-4,ncol(A),ncol(A))), Class="dgCMatrix")
# mix2 <- mmec(Yield~1,
#              random=~vsc(isc(id),Gu=Ai)
#              + Rowf + Colf,
#              rcov=~units,
#              data=DT)
# summary(mix2)$varcomp
#
# vg <- summary(mix2)$varcomp[1,1] # genetic variance
# G <- A*vg # genetic variance-covariance
# Ci <- mix2$Ci # coefficient matrix
# ind <- as.vector(mix2$partitions$`vsc(isc(id), Gu = Ai)`)
# ind <- seq(ind[1],ind[2])
# Ctt <- Ci[ind,ind] # portion of Ci for genotypes
# R2 <- (G - Ctt)/G # reliability matrix
# mean(diag(R2)) # average reliability of the trial
#
# ####====================####
# #### multivariate model ####
# ####     2 traits       ####
# ####====================####
# #### be patient take some time
# ans.m <- mmer(cbind(Yield,color)~1,
#                random=~ vsr(id, Gu=A)
#                + vsr(Rowf,Gtc=diag(2))
#                + vsr(Colf,Gtc=diag(2)),
#                rcov=~ vsr(units),
#                data=DT)
# cov2cor(ans.m$sigma$`u:id`)

}
\keyword{datasets}
\seealso{The core functions of the package \code{\link{mmer}} and \code{\link{mmec}} }

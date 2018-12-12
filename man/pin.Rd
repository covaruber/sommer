\name{pin}
\alias{pin}
\title{pin functionality}
\description{
 Post-analysis procedure to calculate functions of variance components. Its intended use is when the variance components are either simple variances or are variances and covariances in an unstructured matrix. The functions covered are linear combinations of the variance components (for example, phenotypic variance), a ratio of two components (for example, heritabilities) and the correlation based on three components (for example, genetic correlation). 
 
The pin file specifies the functions to be calculated.

The calculations are based on the estimated variance parameters and their variance matrix as represented by the inverse of the Fisher or Average information matrix. Note that this matrix has zero values for fixed variance parameters including those near the parameter space boundary.
}
\usage{
  pin(object, transform)
}

\arguments{
  \item{object}{a model fitted with the mmer or mmer2 functions.}
  \item{transform}{formula to calculate the function.}
}
\value{
  \describe{
    
    \item{$dd}{the parameter and its standard error.}
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
####=========================================####
#### EXAMPLE 1
#### simple example with univariate models
####=========================================####
####=========================================####
# data(CPdata)
# #### create the variance-covariance matrix 
# A <- A.mat(CPgeno)
# #### look at the data and fit the model
# head(CPpheno)
# mix1 <- mmer2(Yield~1,
#               random=~g(id), 
#               G=list(id=A), data=CPpheno)
# summary(mix1)
# #### run the pin function
# pin(mix1, h2 ~ V1 / ( V1 + V2 ) )

# ####=========================================####
# ####=========================================####
# #### EXAMPLE 2
# #### simple example with multivariate models
# ####=========================================####
# ####=========================================####
# data(CPdata)
# CPpheno <- CPdata$pheno
# CPgeno <- CPdata$geno
# #### create the variance-covariance matrix 
# A <- A.mat(CPgeno)
# #### look at the data and fit the model
# head(CPpheno)
# mix2 <- mmer2(cbind(Yield,color)~1,
#                random=~us(trait):g(id), 
#                rcov=~us(trait):units,
#                G=list(id=A), 
#               data=CPpheno, MVM=TRUE)
# summary(mix2)
# ## genetic correlation
# pin(mix2, gen.cor ~ V2 / sqrt(V1*V3))
# 
# ####=========================================####
# ####=========================================####
# #### EXAMPLE 3
# #### more complex multivariate model
# ####=========================================####
# ####=========================================####
# data(BTdata)
# mix3 <- mmer2(cbind(tarsus, back) ~ sex,
#                random = ~ us(trait):dam + diag(trait):fosternest,
#                rcov=~us(trait):units,
#                data = BTdata)
# summary(mix3)
# #### calculate the genetic correlation
# pin(mix3, gen.cor ~ V2 / sqrt(V1*V3))
# 
# ####=========================================####
# ####=========================================####
# #### EXAMPLE 4
# #### going back to simple examples
# ####=========================================####
# ####=========================================####
# data(BTdata)
# mix4 <- mmer2(tarsus ~ sex, random = ~ dam + fosternest,
#                data = BTdata) 
# summary(mix4)
# #### calculate the ratio and its SE
# pin(mix4, dam.prop ~ V1 / ( V1 + V2 + V3 ) )

}
\seealso{The core function of the package \code{\link{mmer}}}
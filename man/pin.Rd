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
# data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix 
# A <- A.mat(GT)
# #### look at the data and fit the model
# head(DT)
# mix1 <- mmer(Yield~1,
#               random=~vs(id,Gu=A), 
#               data=DT)
# summary(mix1)
# #### run the pin function
# pin(mix1, h2 ~ V1 / ( V1 + V2 ) )

# ####=========================================####
# ####=========================================####
# #### EXAMPLE 2
# #### simple example with multivariate models
# ####=========================================####
# ####=========================================####
# data(DT_cpdata)
# DT <- DT_cpdata
# GT <- GT_cpdata
# MP <- MP_cpdata
# #### create the variance-covariance matrix 
# A <- A.mat(GT)
# #### look at the data and fit the model
# head(DT)
# mix2 <- mmer(cbind(Yield,color)~1,
#                random=~vs(id,Gu=A, Gt=unsm(2)), 
#                rcov=~vs(units, Gt=unsm(2)),
#               data=DT)
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
# data(DT_btdata)
# DT <- DT_btdata
# mix3 <- mmer(cbind(tarsus, back) ~ sex,
#                random = ~ vs(dam, Gtc=unsm(2)) + vs(fosternest,Gtc=diag(2)),
#                rcov=~vs(units,Gtc=unsm(2)),
#                data = DT)
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
# data(DT_btdata)
# DT <- DT_btdata
# mix4 <- mmer2(tarsus ~ sex, random = ~ dam + fosternest,
#                data = DT) 
# summary(mix4)
# #### calculate the ratio and its SE
# pin(mix4, dam.prop ~ V1 / ( V1 + V2 + V3 ) )

}
\seealso{The core function of the package \code{\link{mmer}}}
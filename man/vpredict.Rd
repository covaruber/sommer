\name{vpredict}
\alias{vpredict}
\alias{vpredict.mmer}
\title{vpredict form of a LMM fitted with mmer}
\description{\code{vpredict} method for class \code{"mmer"}.

Post-analysis procedure to calculate linear combinations of variance components. Its intended use is when the variance components are either simple variances or are variances and covariances in an unstructured matrix. The functions covered are linear combinations of the variance components (for example, phenotypic variance), a ratio of two components (for example, heritabilities) and the correlation based on three components (for example, genetic correlation). 
 
The calculations are based on the estimated variance parameters and their variance matrix as represented by the inverse of the Fisher or Average information matrix. Note that this matrix has zero values for fixed variance parameters including those near the parameter space boundary.

The transform is specified with a formula. On the left side of the formula is a name for the transformation. On the right side of the formula is a transformation specified with shortcut names like `V1`, `V2`, etc.  The easiest way to identify these shortcut names is to use `summary(object)$varcomp`. The rows of this object can referred to with shortcuts `V1`, `V2`, etc.  See the example below.

}
\usage{

vpredict(object, transform)
\method{vpredict}{mmer}(object, transform)

}
\arguments{

 \item{object}{a model fitted with the mmer function.}
 
 \item{transform}{a formula to calculate the function.}
 
}
\details{
  
  The delta method (e.g., Lynch and Walsh 1998, Appendix 1; Ver Hoef 2012) uses
  a Taylor series expansion to approximate the moments of a function of
  parameters. Here, a second-order Taylor series expansion is implemented to
  approximate the standard error for a function of (co)variance parameters.
  Partial first derivatives of the function are calculated by algorithmic
  differentiation with \code{\link[stats]{deriv}}.

Though \code{vpredict} can calculate standard errors for non-linear functions
  of (co)variance parameters from a fitted \code{mmer} model, it is limited
  to non-linear functions constructed by mathematical operations such as the
  arithmetic operators \code{+}, \code{-}, \code{*}, \code{/} and \code{^},
  and single-variable functions such as  \code{exp} and \code{log}. See 
  \code{\link[stats]{deriv}} for more information.
  
 }
\value{

  \item{dd}{the parameter and its standard error.}
  
 }
 \references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
  Lynch, M. and B. Walsh 1998. Genetics and Analysis of Quantitative Traits. Sinauer Associates, Inc., Sunderland, MA, USA.

  Ver Hoef, J.M. 2012. Who invented the delta method? The American Statistician 66:124-127. DOI: 10.1080/00031305.2012.687494
  
 }
\author{Giovanny Covarrubias}
\seealso{\code{\link{vpredict}}, \code{\link{mmer}}}
\keyword{models}
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
# summary(mix1)$varcomp
# #### run the vpredict function
# vpredict(mix1, h2 ~ V1 / ( V1 + V2 ) )
# 
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
# summary(mix2)$varcomp
# ## genetic correlation
# vpredict(mix2, gen.cor ~ V2 / sqrt(V1*V3))
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
# summary(mix3)$varcomp
# #### calculate the genetic correlation
# vpredict(mix3, gen.cor ~ V2 / sqrt(V1*V3))
# 
# ####=========================================####
# ####=========================================####
# #### EXAMPLE 4
# #### going back to simple examples
# ####=========================================####
# ####=========================================####
# data(DT_btdata)
# DT <- DT_btdata
# mix4 <- mmer(tarsus ~ sex, random = ~ dam + fosternest,
#                data = DT) 
# summary(mix4)$varcomp
# #### calculate the ratio and its SE
# vpredict(mix4, dam.prop ~ V1 / ( V1 + V2 + V3 ) )


}

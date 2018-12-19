\name{bivariateRun}
\alias{bivariateRun}
\title{bivariateRun functionality}
\description{
 Sometimes multi-trait models can present many singularities making the model hard to estimate with many traits. One of the most effective strategies is to estimate all possible variance and covariances splitting in multiple bivariate models. This function takes a model that has t traits and splits the model in as many bivariate models as needed to estimate all the variance and covariances to provide the initial values for the model with all traits.
}
\usage{
  bivariateRun(model, n.core)
}

\arguments{
  \item{model}{a model fitted with the mmer function with argument return.param=TRUE.}
  \item{n.core}{number of cores to use in the mclapply function to parallelize the models to be run to avoid increase in computational time. Please keep in mind that this is only available in Linux and macOS systems. Please check the details in the \code{\link{mclapply}} documentation of the parallel package.}
}
\value{
  \describe{
    
    \item{$sigmas}{the list with the variance covariance parameters for all traits together.}
    \item{$sigmascor}{the list with the correlation for the variance components for all traits together.}
    \item{$model}{the results from the bivariate models.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
  
}
\author{
  Giovanny Covarrubias-Pazaran
}
\examples{

# ####=========================================####
# ####=========================================####
# #### EXAMPLE 1
# #### simple example with univariate models
# ####=========================================####
# ####=========================================####
# data("DT_cpdata")
# #### create the variance-covariance matrix
# A <- A.mat(GT)
# #### look at the data and fit the model
# head(DT)
# ans.m <- mmer(cbind(Yield,color,FruitAver, Firmness)~1,
#                random=~ vs(id, Gu=A, Gtc=unsm(4))
#                + vs(Rowf,Gtc=diag(4))
#                + vs(Colf,Gtc=diag(4)), na.method.Y="include",
#                rcov=~ vs(units,Gtc=unsm(4)), return.param = TRUE,
#                data=DT)
# 
# # define the number of cores (number of bivariate models) as (nt*(nt-1))/2 
# nt=4
# (nt*(nt-1))/2
# res <- bivariateRun(ans.m,n.core = 6)
# # now use the variance componets to fit a join model
# mm <- transformConstraints(ans.m[[8]],3)
# 
# ans.m.final <- mmer(cbind(Yield,color,FruitAver, Firmness)~1,
#               random=~ vs(id, Gu=A, Gtc=unsm(4))
#               + vs(Rowf,Gtc=diag(4))
#               + vs(Colf,Gtc=diag(4)), na.method.Y="include",
#               rcov=~ vs(units,Gtc=unsm(4)), 
#               init = res$sigmas_scaled, constraints = mm,
#               data=DT, iters=1)
# 
# summary(ans.m.final)

}
\seealso{The core function of the package \code{\link{mmer}}}
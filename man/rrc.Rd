\name{rrc}
\alias{rrc}
\title{reduced rank covariance structure}
\description{
  \code{rrc} creates a reduced rank factor analytic covariance structure selecting the n vectors of the L matrix of the Cholesky decomposition or the U vectors of the SVD decomposition to create a new incidence matrix that can be used with the \code{\link{mmec}} solver.
}
\usage{
  rrc(timevar=NULL, idvar=NULL, response=NULL, 
      Gu=NULL, nPC=2, returnGamma=FALSE, cholD=TRUE)
}

\arguments{
  \item{timevar}{vector of the dataset containing the variable to be used to form columns in the wide table.}
  \item{idvar}{vector of the dataset containing the variable to be used to form rows in the wide table.}
  \item{response}{vector of the dataset containing the response variable to be used to fill the cells of the wide table.}
  \item{Gu}{an optional covariance matrix (\strong{NOT THE INVERSE}) between levels of the idvar in case a sparse (unbalanced) design between timevar and idvar exist.}
  \item{nPC}{number of principal components to keep.}
  \item{returnGamma}{a TRUE/FALSE argument specifying if the function should return the matrix of loadings of the incidence matrix for the model. The default is FALSE so it returns the incidence matrix.}
    \item{cholD}{a TRUE/FALSE argument specifying if the Cholesky decomposition should be calculated or the singular value decomposition should be used instead.}
}
\details{
This implementation of a version of the reduced rank factor analytic models uses the so-called principal component (PC) models (Meyer, 2009) which assumes specific effects (psi) are equal to 0. The model is as follows:

y = Xb + Zu + e

where the variance of u ~ MVN(0, Sigma)

Sigma = (Gamma_t Gamma) + Psi

\strong{Extended factor analytic model}:

y = Xb + Z(I Gamma)c + Zs + e = Xb + Z*c + Zs + e

where \code{y} is the response variable, \code{X} and \code{Z} are incidence matrices for fixed and random effects respectively, \code{I} is a diagonal matrix, \code{Gamma} are the factor loadings for \code{c} common factor scores, and \code{s} are the specific effects, \code{e} is the vector of residuals.   

\strong{Reduced rank model:}

y = Xb + Z(I Gamma)c + e = Xb + Z*c + e

which is equal to the one above but assumes specific effects = 0. 

\strong{The algorithm is the following:}

1) creates a wide-format table of timevar (m columns) by idvar (q rows)  named H that is used to form the initial covariance matrix Sigma = H'H of dimensions m x m. 

2) The Sigma matrix is then center and scaled. 

3) A Cholesky (L matrix) or SVD decomposition (U D V') is performed in the Sigma matrix.

4) n vectors from L (when Cholesky is used) or U sqrt(D) (when SVD is used) are kept to form Gamma.

4) Gamma is used to form a new incidence matrix as Z* = Z Gamma

5) This matrix is later used for the REML machinery to be used with the usc (unstructured) or dsc (diagonal) structures to estimate variance components and factor scores. The resulting BLUPs from the mixed model are the optimized factor scores.

This implementation does not update the loadings during the REML process, only estimates the REML factor scores for fixed loadings. This is different to other software (e.g., asreml) where the loadings are updated as well. 

BLUPs for genotypes in all locations can be recovered as:

u = Gamma u_scores

The resulting loadings (Gamma) and factor scores can be thought as an equivalent to the classical factor analysis.

}
\value{
  \describe{
    
    \item{$Z}{ a incidence matrix Z* = Z Gamma which is the original incidence matrix for the timevar multiplied by the loadings.}
  }
}
\references{

  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744
  
  Meyer K (2009) Factor analytic models for genotype by environment type problems and structured covariance matrices. Genetics Selection Evolution, 41:21
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_h2)
# DT <- DT_h2
# DT=DT[with(DT, order(Env)), ]
# 
# ans1b <- mmec(y~Env,
#               random=~vsc( usc(rrc(Env, Name, y, nPC = 3)) , isc(Name)),
#               rcov=~units, 
#               # we recommend giving more iterations to these models
#               nIters = 50, 
#               # we recommend giving more EM iterations at the beggining 
#               emWeight = c(rep(1,10),logspace(10,1,.05), rep(.05,80)),
#               data=DT)
# 
# summary(ans1b)$varcomp
# vcd <- diag(ans1b$theta[[1]])
# vcd/sum(vcd) # the 3rd PC still explains more than the 2nd
# ## Extract BLUPs
# ## extract loadings
# loadings=with(DT, rrc(Env, Name, y, returnGamma = TRUE, nPC = 3))$Gamma
# ## extract factor scores
# scores <- ans1b$uList[[1]]; 
# ## BLUPs for all environments
# E= scores  t(loadings)  # change space for matrix product
# ## Extract the covariance matrix
# vc <- ans1b$theta[[1]];vc
# G = loadings  vc  t(loadings) # change space for matrix product
# lattice::levelplot(cov2cor(G))


}
\seealso{The function \code{\link{vsc}} to know how to use \code{rrc} in the \code{\link{mmec}} solver.}
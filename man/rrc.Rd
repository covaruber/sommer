\name{rrc}
\alias{rrc}
\title{reduced rank covariance structure}
\description{
  \code{rrc} creates a reduced rank factor analytic covariance structure selecting the n principal components (nPC) to be used with the \code{\link{mmec}} solver.
}
\usage{
  rrc(timevar=NULL, idvar=NULL, response=NULL, 
      Gu=NULL, nPC=2, returnLam=FALSE, cholD=TRUE)
}

\arguments{
  \item{timevar}{vector of the dataset containing the variable to be used to form columns in the wide table.}
  \item{idvar}{vector of the dataset containing the variable to be used to form rows in the wide table.}
  \item{response}{vector of the dataset containing the response variable to be used to fill the cells of the wide table.}
  \item{Gu}{an optional covariance matrix (\strong{NOT THE INVERSE}) between levels of the idvar in case a sparse (unbalanced) design between timevar and idvar exist.}
  \item{nPC}{number of principal components to keep.}
  \item{returnLam}{a TRUE/FALSE argument specifying if the function should return the matrix of loadings of the incidence matrix for the model. The default is FALSE so it returns the incidence matrix.}
  \item{cholD}{a TRUE/FALSE argument specifying if the function should perform a Cholesky decomposition (TRUE) to form the matrix of loadings or the singular value decomposition (FALSE). Default is Cholesky.}
}
\details{
This implementation of a version of the reduced rank factor analytic models uses the so-called principal component (PC) models (Meyer, 2009) which assumes specific effects (psi) are equal to 0. The function does the following:

1) creates a wide-format table of timevar (columns) by idvar (rows) to form the initial covariance matrix Sigma. 

2) The Sigma matrix is then center and scaled. 

3) A Cholesky or SVD decomposition is performed in the Sigma matrix to return the matrix of loadings Lam.

4) n PCs are kept and used to form a new incidence matrix Z.star = Z Lam

5) This matrix is later used for the mmec() machinery to be used with the usc (FA) structure. The resulting BLUPs from the mixed model are the optimized factor scores.

This implementation does not update the loadings during the REML process, only estimates the factor scores. This is different to the asreml FA modeling where the loadings are estimated as well. 

}
\value{
  \describe{
    
    \item{$Z}{ a incidence matrix Z.star = Z Lam which is the original incidence matrix for the timevar multiplied by the loadings.}
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
DT <- DT_h2
DT=DT[with(DT, order(Env)), ]

ans1b <- mmec(y~Env,
              random=~vsc( usc(rrc(Env, Name, y, nPC = 3)) , isc(Name)),
              rcov=~units, 
              # we recommend giving more iterations to these models
              nIters = 50, 
              # we recommend giving more EM iterations at the beggining 
              emWeight = c(rep(1,10),logspace(10,1,.05), rep(.05,80)),
              data=DT)

summary(ans1b)$varcomp
vcd <- diag(ans1b$theta[[1]])
vcd/sum(vcd) # the 3rd PC still explains more than the 2nd
# ## Extract BLUPs
# ## extract loadings
# loadings=with(DT, rrc(Env, Name, y, returnLam = TRUE, nPC = 3))$Lam 
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
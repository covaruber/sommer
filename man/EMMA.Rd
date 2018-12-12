\name{MEMMA}
\alias{MEMMA}
\title{Multivariate Efficient Mixed Model Association Algorithm}
\description{
This function is used internally in the function \code{\link{mmer}} when multiple responses are selected for a single variance component other than the error. It uses the efficient mixed model association (MEMMA) algorithm. 
}
\usage{
MEMMA(Y, X=NULL, ZETA=NULL, tolpar = 1e-06, tolparinv = 1e-06, check.model=TRUE,
      silent=TRUE)
}
\arguments{
  \item{Y}{a numeric vector for the response variable}
  \item{X}{an incidence matrix for fixed effects.}
  \item{ZETA}{an incidence matrix for random effects. This can be for one or more random effects. This NEEDS TO BE PROVIDED AS A LIST STRUCTURE. For example Z=list(list(Z=Z1, K=K1), list(Z=Z2, K=K2), list(Z=Z3, K=K3)) makes a 2 level list for 3 random effects. The general idea is that each random effect with or without its variance-covariance structure is a list, i.e. list(Z=Z1, K=K1) where Z is the incidence matrix and K the var-cov matrix. When moving to more than one random effect we need to make several lists that need to be inside another list. What we call a 2-level list, i.e. list(Z=Z1, K=K1) and list(Z=Z2, K=K2) would need to be put in the form; list(list(Z=Z1, K=K1),list(Z=Z1, K=K1)), which as can be seen, is a list of lists (2-level list). }
  \item{tolpar}{tolerance parameter for convergence}
  \item{tolparinv}{tolerance parameter for matrix inverse}
  \item{check.model}{a TRUE/FALSE value indicating if list structure provided by the user is correct to fix it. The default is TRUE but is turned off to FALSE within the mmer function which would imply a double check.}
  \item{silent}{a TRUE/FALSE value indicating if the function should draw the progress bar or iterations performed while working or should not be displayed.}
  
}
\details{

.

The likelihood function optimized in this algorithm is:

.

logL = (n - p) * log(sum(eta^2/{ lambda + delta})) + sum(log(lambda + delta))

.

where: 
(n-p)  refers to the degrees of freedom
lambda are the eigenvalues mentioned by Kang et al.(2008) 
delta is the REML estimator of the ridge parameter

.

The algorithm can be summarized in the next steps:

.

1) provide initial value for the ridge parameter

2) estimate S = I - X(X'X)-X'

3) obtain the phenotypic variance V = ZKZ' + delta.prov*I 

4) perform an eigen decomposition of SVS 

5) create "lambda"" as the eigenvalues of SVS and "U"" as the eigenvectors

6) estimate eta=U'y 

7) optimize the likelihood shown above providing "eta", "lambdas" and optimize with respect to "delta" which is the ridge parameter and contains Ve/Vu


}
\value{
If all parameters are correctly indicated the program will return a list with the following information:
\describe{

\item{$Vu}{ a scalar value for the variance component estimated}
\item{$Ve}{ a scalar value for the error variance estimated}
\item{$V.inv}{ a matrix with the inverse of the phenotypic variance V = ZGZ+R, V^-1}
\item{$u.hat}{ a vector with BLUPs for random effects}
\item{$Var.u.hat}{ a vector with variances for BLUPs}
\item{$PEV.u.hat}{ a vector with predicted error variance for BLUPs}
\item{$beta.hat}{ a vector for BLUEs of fixed effects}
\item{$Var.beta.hat}{ a vector with variances for BLUEs}
\item{$X}{ incidence matrix for fixed effects, if not passed is assumed to only include the intercept}
\item{$Z}{ incidence matrix for random effects, if not passed is assumed to be a diagonal matrix}
\item{$K}{ the var-cov matrix for the random effect fitted in Z}
\item{$ll}{ the log-likelihood value for obtained when optimizing the likelihood function when using ML or REML}

}
}
\references{

Kang et al. 2008. Efficient control of population structure in model organism association mapping. Genetics 178:1709-1723.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 

}
\examples{

####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####
# data(CPdata)
# ### look at the data
# head(CPpheno)
# CPgeno[1:5,1:5]
# ## fit a model including additive and dominance effects
# Y <- CPpheno[,c("color","Yield")]
# Za <- diag(dim(Y)[1])
# A <- A.mat(CPgeno) # additive relationship matrix
# ####================####
# #### ADDITIVE MODEL ####
# ####================####
# ETA.A <- list(add=list(Z=Za,K=A))
# #ans.A <- MEMMA(Y=Y, ZETA=ETA.A)
# #ans.A$var.comp
}
\seealso{The core functions of the package \code{\link{mmer}} }

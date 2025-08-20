\name{MNR}
\alias{MNR}

\title{
Multivariate Newton-Raphson algorithm 
}
\description{
Multivariate Newton-Raphson algorithm used behind mmes when the henderson argument 
is set to FALSE. Algorithm made available for users that want to avoid the user-friendly interface 
and provide their matrices directly.
}
\usage{
MNR(Y, 
    X, Gx,
    Z, K,
    R, 
    Ge, GeI, 
    W, isInvW,
    iters, tolpar, tolparinv,
    ai, pev,
    verbose, retscaled,
    stepweight, emweight
    )
}

\arguments{
  \item{Y}{
A matrix with rows for records and columns for traits. Expected class is a 'matrix'.
}
  \item{X}{
Design matrices for fixed effects. Expected class is a list with as many design 
matrices as fixed effects.
}
  \item{Gx}{
Constraints for fixed effects. Expected class is a list with as many matrices as 
fixed effects. Each matrix is a square matrix with as many rows and columns as
levels of the fixed effects.
}
  \item{Z}{
Design matrices for random effects. Expected class is a list with as many design 
matrices as random effects. Each element in the list should be a matrix of class
'dgCMatrix'
}
  \item{K}{
Covariance matrices for design matrices of random effects. Expected class is a list with as many  
covariance matrices as random effects specified in Z. Each element in the list should be a matrix of class
'dgCMatrix'
}
  \item{R}{
Residual matrices for residual effects. Expected class is a list with as many  
residual matrices as residual effects. Each element in the list should be a square 
matrix of dimensions n x n, where n is the number of records. Each matrix should be of class
'matrix'
}
  \item{Ge}{
Initial values for variance components. Expected class is a list with as many  
variance component matrices as random effects specified in Z. Each element in the 
list should be a matrix of dimensions t x t, where t is the number of traits and of 
class 'matrix'. Initial values are any real values.
}
  \item{GeI}{
Initial constraints for variance components. Expected class is a list with as many  
variance component constraints matrices as random effects specified in Z. Each element in the 
list should be a matrix of dimensions t x t, where t is the number of traits and of 
class 'matrix'. Values expected are:

0: not to be estimated

1: estimated and constrained to be positive (i.e. variance component)

2: estimated and unconstrained (can be negative or positive, i.e. covariance component)

3: not to be estimated but fixed (value has to be provided in the Gti argument)

}
  \item{W}{
Design matrix for weights. A matrix of class 'matrix' for weighting the records.
}
  \item{isInvW}{
A value of class 'logical' to indicate if the W matrix provided is already an inverse or
not. This aims to speed up the computations.
}
  \item{iters}{
Maximum number of iterations allowed in REML. Value is of class 'integer'.
}
  \item{tolpar}{
Convergence criteria for the change in log-likelihood. Value is of class 'numeric'.
}
  \item{tolparinv}{
Tolerance parameter for matrix inverse used when singularities are encountered in
the estimation procedure. Value is of class 'numeric'.
}
  \item{ai}{
A value of class 'logical' to indicate if the Average Information algorithm
should be used instead. That is faster but much less stable.
}
  \item{pev}{
A value of class 'logical' to indicate if the predicted error variance should be 
computed or not. If FALSE computations are speeded up for models with many effects
to be estimated.
}
  \item{verbose}{
A value of class 'logical' to value to indicate if the program should return
the progress of the iterative algorithm.
}
  \item{retscaled}{
A value of class 'logical' to indicate if we should avoid scaling the traits.
}
  \item{stepweight}{
Vector of class 'numeric' with length equal to niters to specify the relative weight
given to the second derivative Newton update.
}
  \item{emweight}{
Vector of class 'numeric' with length equal to niters to specify the relative weight
given to the EM update compared to the NR update.
}


}
\details{

This is the Rcpp-coded Direct-Inversion LMM REML algorithm used when the 'henderson'
argument is set to FALSE.

}
\value{
If all parameters are correctly indicated the program will return a list with the following information:

\item{Vi}{the inverse of the phenotypic variance matrix V^- = (ZGZ+R)^-1}
\item{P}{the projection matrix  Vi - [Vi*(X*Vi*X)^-*Vi] }
\item{sigma}{a list with the values of the variance-covariance components with one list element for each random effect.}
\item{sigma_scaled}{a list with the values of the scaled variance-covariance components with one list element for each random effect.}
\item{sigmaSE}{Hessian matrix containing the variance-covariance for the variance components. SE's can be obtained taking the square root of the diagonal values of the Hessian.}
\item{Beta}{a data frame for trait BLUEs (fixed effects).}
\item{VarBeta}{a variance-covariance matrix for trait BLUEs}
\item{U}{a list (one element for each random effect) with a data frame for trait BLUPs.}
\item{VarU}{a list (one element for each random effect) with the variance-covariance matrix for trait BLUPs.}
\item{PevU}{a list (one element for each random effect) with the predicted error variance matrix for trait BLUPs.}
\item{fitted}{ Fitted values y.hat=XB}
\item{residuals}{Residual values e = Y - XB}
\item{AIC}{Akaike information criterion}
\item{BIC}{Bayesian information criterion}
\item{convergence}{a TRUE/FALSE statement indicating if the model converged.}
\item{monitor}{The values of log-likelihood and variance-covariance components across iterations during the REML estimation.}
\item{percChange}{The percent change of variance components across iterations. There should be one column less than the number of iterations. Calculated as percChange = ((x_i/x_i-1) - 1) * 100 where i is the ith iteration.}
\item{dL}{The vector of first derivatives of the likelihood with respect to the ith variance-covariance component.}
\item{dL2}{The matrix of second derivatives of the likelihood with respect to the i.j th variance-covariance component.}

}

\references{

Mrode, R. A. (2014). Linear models for the prediction of animal breeding values. Cabi.

Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744

}

\examples{

data(DT_cpdata)
DT <- DT_cpdata

# response matrix
y <- matrix(imputev(DT$Yield))
# fixed effect incidence matrix
X <- model.matrix(~Rowf, data=DT)
# random effect incidence matrix
Z <- Matrix::sparse.model.matrix(~id - 1, data=DT)
colnames(Z) <- gsub("id","",colnames(Z))
# covariance for random effect
GT <- GT_cpdata
A <- A.mat(GT) # additive relationship matrix
A <- A[colnames(Z),colnames(Z)] # make sure of the order
# residual effect incidence matrix
R1 <- Matrix::Diagonal(n=length(y))
# weights matrix
W <- diag(length(y))
# some other parameters
maxIter=10
stepWeight <- rep(0.9, maxIter)
stepWeight[1:2] <- c(0.5, 0.7)
emWeights <- rep(0,maxIter)
# model fit
res <- MNR( Y=matrix(y), # response
                     X=list(X), Gx=list(matrix(1)), # fixed effects
                     Z=list(Z), K=list(A), # random effects
                     R=list(R1), # residual effects
                     Ge=list(matrix(0.15),matrix(0.75)), # inital vc 
                     GeI=list(matrix(1),matrix(1)), # vc constraints
                     W=W, isInvW=TRUE, # weights for records
                     iters=maxIter,
                     tolpar=1e-4, tolparinv=1e-6, # tolerance
                     ai=FALSE, pev=FALSE, # algorithm specifics
                     verbose=TRUE, 
                     retscaled=FALSE, 
                     stepweight=stepWeight, # second derivatives weights
                     emweight=emWeights # em update weights
)
res$sigma
res$Beta
res$U[[1]]

}
\seealso{\code{\link{mmes}} -- the core function of the package }

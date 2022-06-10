\name{AI}
\alias{AI}
\title{Average Information Algorithm}
\description{
Univariate version of the average information (AI) algorithm.
}
\usage{
AI(X=NULL,Z=NULL, Zind=NULL, Ai=NULL,y=NULL,S=NULL, H=NULL,
    nIters=80, tolParConvLL=1e-4, tolParConvNorm=.05, tolParInv=1e-6, theta=NULL, 
    thetaC=NULL, thetaF=NULL, addScaleParam=NULL,
    weightInfEMv=NULL, weightInfMat=NULL)
}
\arguments{
    \item{X}{an incidence matrix for fixed effects.}
  \item{Z}{Z is a list of lists each element contains the Z matrices required for the covariance structure specified for a random effect.}
  \item{Zind}{vector specifying to which random effect each Z matrix belongs to. }
  \item{Ai}{is a list with the inverses of the relationship matrix for each random effect.}
  \item{y}{is the response variable}
  \item{S}{is the list of residual matrices.}
  \item{H}{is the matrix of weights. This will be squared via the cholesky decomposition and apply to the residual matrices.}
  \item{nIters}{number of REML iterations .}
  \item{tolParConvLL}{rule for stoping the optimization problem, difference in log-likelihood between the current and past iteration.}
  \item{tolParConvNorm}{rule for stoping the optimization problem, difference in norms.}
  \item{tolParInv}{value to add to the diagonals of a matrix that cannot be inverted because is not positive-definite.}
  \item{theta}{ is the initial values for the vc (matrices should be symmetric).}
  \item{thetaC}{ is the constraints for vc: 1 positive, 2 unconstrained, 3 fixed.}
  \item{thetaF}{ is the dataframe indicating the fixed constraints as x times another vc, rows indicate the variance components, columns the scale parameters (other VC plus additional ones preferred).}
  \item{addScaleParam}{ any additional scale parameter to be included when applying constraints in thetaF.}
  \item{weightInfEMv}{ is the vector to be put in a diagonal matrix (a list with as many matrices as iterations) representing the weight assigned to the EM information matrix.}
  \item{weightInfMat}{ is a vector of weights to the information matrix for the operation delta = I- * dLu/dLx # unstructured models may require less weight to the information matrix.}
}
\details{
This algorithm is based on Jensen, Madsen and Thompson (1997) 
}
\value{
If all parameters are correctly indicated the program will return a list with the following information:
\describe{
\item{res}{ a list of different outputs}
}
}
\references{
Jensen, J., Mantysaari, E. A., Madsen, P., and Thompson, R. (1997). Residual maximum likelihood estimation of (co) variance components in multivariate mixed linear models using average information. Journal of the Indian Society of Agricultural Statistics, 49, 215-236.

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744   
}
\examples{
####=========================================####
#### For CRAN time limitations most lines in the 
#### examples are silenced with one '#' mark, 
#### remove them and run the examples
####=========================================####

data("DT_example")
DT <- DT_example
K <- A_example
#### look at the data and fit the model
head(DT)

zz <- with(DT, vsr(dsr(Env),Name))

Z <- c(list(model.matrix(~Name-1, data=DT)),zz$Z)

Zind <- c(1,2,2,2)

A <- list(diag(41), diag(41))#rep(list(diag(41)),4)

Ai <- lapply(A, function(x){solve(x)})

theta <- list(
  matrix(10,1,1),
  diag(10,3,3),
  diag(10,3,3)
);theta

thetaC <- list(
  matrix(1,1,1),
  diag(1,3,3),
  diag(1,3,3)
);thetaC

X <- model.matrix(~Env, data=DT)

y <- as.matrix(DT$Yield)

DTx <- DT; DTx$units <- as.factor(1:nrow(DTx))
ss <- with(DTx, vsr(dsr(Env),units) )

S <- ss$Z 

H <- diag(length(y))

addScaleParam <- 0
nn <- unlist(lapply(thetaC, function(x){length(which(x > 0))}))
nn2 <- sum(nn[1:max(Zind)])
ff <- diag(nn2)
thetaF <- cbind(ff,matrix(0,nn2,1))

## apply the function
weightInfMat=rep(1,40); # weights for the information matrix
weightInfEMv=c(seq(.9,.1,-.1),rep(0,36)); # weights for the EM information matrix

# expr = res3<-AI(X=X,Z=Z, Zind=Zind,
#                 Ai=Ai,y=y,
#                 S=S,
#                 H=H, 
#                 nIters=20, tolParConvLL=1e-5,
#                 tolParConvNorm=0.05,
#                 tolParInv=1e-6,theta=theta,
#                 thetaC=thetaC,thetaF=thetaF,
#                 addScaleParam=addScaleParam, weightInfEMv = weightInfEMv,
#                 weightInfMat = weightInfMat
#                 
# )
# # compare results
# res3$monitor


}
\seealso{The core functions of the package \code{\link{mmer}} }
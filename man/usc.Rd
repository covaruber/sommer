\name{usc}
\alias{usc}
\title{unstructured covariance structure}
\description{
  \code{usc} creates an unstructured covariance structure for specific levels of the random effect to be used with the \code{\link{mmec}} solver.
}
\usage{
  usc(x, thetaC, theta)
}

\arguments{
  \item{x}{vector of observations for the random effect.}
    \item{thetaC}{an optional symmetric matrix for constraints in the variance-covariance components. The symmetric matrix should have as many rows and columns as the number of levels in the factor 'x'. The values in the matrix define how the variance-covariance components should be estimated:
  
  0: component will not be estimated
  
  1: component will be estimated and constrained to be positive
  
  2: component will be estimated and unconstrained
  
  3: component will be fixed to the value provided in the theta argument
  
  }
  \item{theta}{an optional symmetric matrix for initial values of the variance-covariance components. The symmetric matrix should have as many rows and columns as the number of levels in the factor 'x'. The values in the matrix define the initial values of the variance-covariance components that will be subject to the constraints provided in thetaC. If not provided, initial values will be calculated as:
  
  diag(ncol(mm))*.05 + matrix(.1,ncol(mm),ncol(mm))
  
  where mm is the incidence matrix for the factor 'x'.
  
  }
}
\value{
  \describe{
    
    \item{$res}{ a list with the provided vector and the variance covariance structure expected for the levels of the random effect.}
  }
}
\references{
  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744 
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{
x <- as.factor(c(1:5,1:5,1:5));x
usc(x)
## how to use the theta and thetaC arguments:
# data(DT_example)
# DT <- DT_example
# theta <- matrix(9:1,3,3); 
# theta[lower.tri(theta)] <- t(theta)[lower.tri(theta)] 
# theta # initial VCs
# thetaC <- fixm(3); thetaC # fixed VCs
# ans1 <- mmec(Yield~Env,
#              random= ~ vsc( usc(Env,theta = theta,thetaC = thetaC),isc(Name) ),
#              rcov= ~ units, nIters = 1,
#              data=DT)
# summary(ans1)$varcomp

}
\seealso{The function \code{\link{vsc}} to know how to use \code{usc} in the \code{\link{mmec}} solver.}
\name{atm}
\alias{atm}
\title{atm covariance structure}
\description{
  \code{atm} creates a diagonal covariance structure for specific levels of the random effect to be used with the \code{\link{mmes}} solver.
}
\usage{
  atm(x, levs, thetaC, theta)
}

\arguments{
  \item{x}{vector of observations for the random effect.}
  \item{levs}{levels of the random effect to use for building the incidence matrices.}
  \item{thetaC}{an optional symmetric matrix for constraints in the variance-covariance components. The symmetric matrix should have as many rows and columns as the number of levels in the factor 'x'. The values in the matrix define how the variance-covariance components should be estimated:
  
  0: component will not be estimated
  
  1: component will be estimated and constrained to be positive
  
  2: component will be estimated and unconstrained
  
  3: component will be fixed to the value provided in the theta argument
  
  }
  \item{theta}{an optional symmetric matrix for initial values of the variance-covariance components. The symmetric matrix should have as many rows and columns as the number of levels in the factor 'x'. The values in the matrix define the initial values of the variance-covariance components that will be subject to the constraints provided in thetaC. If not provided, initial values will be calculated as:
  
  theta* = diag(ncol(mm))*.05 + matrix(.1,ncol(mm),ncol(mm))
  
  where mm is the incidence matrix for the factor 'x'. The values provided should be scaled by the variance of the response variable.
  
  theta = theta*/var(y)
  
  }
}
\value{
  \describe{
    
    \item{$res}{ a list with the provided vector and the variance covariance structure expected.}
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
atm(x, c("1","2"))
## how to use the theta and thetaC arguments:
# data(DT_example)
# DT <- DT_example
# theta <- diag(2)*2; theta # initial VCs
# thetaC <- diag(2)*3; thetaC # fixed VCs
# ans1 <- mmes(Yield~Env,
#              random= ~ vsm( atm(Env, levs=c("CA.2013", "CA.2011"),
#                             theta = theta,thetaC = thetaC),ism(Name) ),
#              rcov= ~ units, nIters = 1,
#              data=DT)
# summary(ans1)$varcomp

}

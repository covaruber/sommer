\name{covm}
\alias{covm}
\title{covariance between random effects}
\description{
  \code{covm} merges the incidence matrices and covariance matrices of two random effects to fit an unstructured model between 2 different random effects to be fitted with the \code{\link{mmes}} solver.
}
\usage{
  covm(ran1, ran2, thetaC=NULL, theta=NULL)
}

\arguments{
  \item{ran1}{the random call of the first random effect.}
  \item{ran2}{the random call of the first random effect.}
  \item{thetaC}{an optional matrix for constraints in the variance components.}
    \item{theta}{an optional symmetric matrix for initial values of the variance-covariance components. When providing customized values, these values should be scaled with respect to the original variance. For example, to provide an initial value of 1 to a given variance component, theta would be built as:
  
  theta = matrix( 1 / var(response) ) 
  
  The symmetric matrix should have as many rows and columns as the number of levels in the factor 'x'. The values in the matrix define the initial values of the variance-covariance components that will be subject to the constraints provided in thetaC. If not provided, initial values will be calculated as:
  
  theta = diag(ncol(mm))*.05 + matrix(.1,ncol(mm),ncol(mm))
  
    where mm is the incidence matrix for the factor 'x'.
  
  }
}
\details{
This implementation aims to fit models where covariance between random variables is expected to exist. For example, indirect genetic effects.

}
\value{
  \describe{
    
    \item{$Z}{ a incidence matrix Z* = Z Gamma which is the original incidence matrix for the timevar multiplied by the loadings.}
  }
}
\references{

  Covarrubias-Pazaran G (2016) Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE 11(6): doi:10.1371/journal.pone.0156744
  
  Bijma, P. (2014). The quantitative genetics of indirect genetic effects: a selective review of modelling issues. Heredity, 112(1), 61-69.
  
}
\author{
  Giovanny Covarrubias-Pazaran
}

\examples{

data(DT_ige)
DT <- DT_ige
covRes <- with(DT, covm( vsm(ism(focal)) , vsm(ism(neighbour)) ) )
str(covRes)
# look at DT_ige help page to see how to fit an actual model

}
\seealso{The function \code{\link{vsm}} to know how to use \code{covm} in the \code{\link{mmes}} solver.}

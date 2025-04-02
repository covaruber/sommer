\name{ism}
\alias{ism}
\title{identity covariance structure}
\description{
  \code{ism} creates an identity covariance structure for the levels of the random effect  to be used with the \code{\link{mmes}} solver. Any random effect with a special covariance structure should end with an ism() structure.
}
\usage{
  ism(x, thetaC=NULL, theta=NULL)
}

\arguments{
  \item{x}{vector of observations for the random effect.}
    \item{thetaC}{an optional 1 x 1 matrix for constraints in the variance-covariance components. The values in the matrix define how the variance-covariance components should be estimated:
  
  0: component will not be estimated
  
  1: component will be estimated and constrained to be positive (default)
  
  2: component will be estimated and unconstrained
  
  3: component will be fixed to the value provided in the theta argument
  
  }
  \item{theta}{an optional 1 x 1 matrix for initial values of the variance-covariance component. When providing customized values, these values should be scaled with respect to the original variance. For example, to provide an initial value of 1 to a given variance component, theta would be built as:
  
  theta = matrix( 1 / var(response) )
  
  The values in the matrix define the initial values of the variance-covariance components that will be subject to the constraints provided in thetaC. If not provided, initial values (theta) will be 0.15
  
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
ism(x)

# data(DT_example)
# ans1 <- mmes(Yield~Env,
#              random= ~ vsm( ism( Name ) ),
#              data=DT_example)
# summary(ans1)$varcomp

}
\seealso{See the function \code{\link{vsm}} to know how to use \code{ism} in the \code{\link{mmes}} solver.}
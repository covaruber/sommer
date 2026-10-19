\name{usm}
\alias{usm}
\title{Unstructured positive-definite covariance structure}
\description{
Unstructured positive-definite covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
usm(x, theta = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Variable or design defining q covariance levels.}
  \item{theta}{Optional positive-definite q by q starting covariance matrix.}
  \item{fixed}{Logical vector of length q(q+1)/2 - 1 controlling the estimable normalized-Cholesky parameters.}
}
\details{
The covariance shape is parameterized as \deqn{K=LL^{\mathsf T},} with lower-triangular \eqn{L} and \eqn{L_{11}=1}.  Diagonal elements after the first are positive and represented on log scales; lower off-diagonal elements are unrestricted.  This guarantees positive definiteness without repeatedly repairing an unconstrained covariance matrix.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(usm(environment), ism(genotype))
}
}
\seealso{
\code{\link{vsm}}, \code{\link{mmes}}.
}

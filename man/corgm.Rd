\name{corgm}
\alias{corgm}
\title{General positive-definite correlation structure}
\description{
General positive-definite correlation structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
corgm(x, theta = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{theta}{Optional positive-definite q by q starting covariance or correlation matrix.}
  \item{fixed}{Logical vector of length q(q-1)/2.}
}
\details{
An unrestricted SPD correlation matrix is represented using a unit-diagonal lower factor \eqn{A}.  The intermediate matrix \eqn{S=AA^{\mathsf T}} is standardized to correlation scale, \deqn{K_{ij}=S_{ij}/\sqrt{S_{ii}S_{jj}}.} The unrestricted working parameters are the strict lower-triangular entries of \eqn{A}.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(corgm(environment), ism(genotype))
}
}
\seealso{
\code{\link{corhm}}, \code{\link{usm}}, \code{\link{vsm}}.
}

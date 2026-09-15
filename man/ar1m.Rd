\name{ar1m}
\alias{ar1m}
\alias{arm}
\title{First-order autoregressive covariance structure}
\description{
First-order autoregressive covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
ar1m(x, rho = 0.30, fixed = FALSE)
arm(x, rho = 0.30, fixed = FALSE)
}
\arguments{
  \item{x}{Ordered factor or design defining q ordered covariance levels.}
  \item{rho}{Starting AR(1) correlation in (-1,1).}
  \item{fixed}{Logical indicating whether rho is fixed.}
}
\details{
The AR(1) shape is \deqn{K_{ij}=\rho^{|i-j|}.} The working coordinate is \eqn{\eta=\operatorname{atanh}(\rho)}, ensuring \eqn{|\rho|<1}.  \code{arm} is an alias for \code{ar1m}.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(ar1m(time), ism(id))
}
}
\seealso{
\code{\link{ar2m}}, \code{\link{ar3m}}, \code{\link{toeplitzm}}, \code{\link{vsm}}.
}

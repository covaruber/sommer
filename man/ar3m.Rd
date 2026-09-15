\name{ar3m}
\alias{ar3m}
\title{Third-order autoregressive covariance structure}
\description{
Third-order autoregressive covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
ar3m(x, pacf = c(0.20, 0.10, 0.05), fixed = NULL)
}
\arguments{
  \item{x}{Ordered factor defining covariance levels.}
  \item{pacf}{Length-three vector of starting partial autocorrelations, each strictly between -1 and 1.}
  \item{fixed}{Logical vector of length three.}
}
\details{
A stationary AR(3) covariance is parameterized through three partial autocorrelations.  The PACF parameterization guarantees stationarity while allowing unconstrained working coordinates.  Reported parameters are the PACFs rather than ordinary AR coefficients.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(ar3m(time), ism(id))
}
}
\seealso{
\code{\link{ar1m}}, \code{\link{ar2m}}, \code{\link{toeplitzm}}, \code{\link{vsm}}.
}

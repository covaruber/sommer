\name{ar2m}
\alias{ar2m}
\title{Second-order autoregressive covariance structure}
\description{
Second-order autoregressive covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
ar2m(x, pacf = c(0.20, 0.10), fixed = NULL)
}
\arguments{
  \item{x}{Ordered factor defining covariance levels.}
  \item{pacf}{Length-two vector of starting partial autocorrelations, each strictly between -1 and 1.}
  \item{fixed}{Logical vector of length two.}
}
\details{
A stationary AR(2) covariance is parameterized through reflection coefficients / partial autocorrelations.  Each PACF is mapped from an unconstrained working coordinate by \code{tanh}, then converted to stable AR coefficients by the Levinson/Durbin recursion.  The PACFs are the reported natural parameters.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(ar2m(time), ism(id))
}
}
\seealso{
\code{\link{ar1m}}, \code{\link{ar3m}}, \code{\link{toeplitzm}}, \code{\link{vsm}}.
}

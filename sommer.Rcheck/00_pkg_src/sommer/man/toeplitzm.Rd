\name{toeplitzm}
\alias{toeplitzm}
\title{General positive-definite Toeplitz correlation structure}
\description{
General positive-definite Toeplitz correlation structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
toeplitzm(x, pacf = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Ordered factor or design defining q covariance levels.}
  \item{pacf}{Optional q-1 starting reflection coefficients / partial autocorrelations.  Defaults to 0.10.}
  \item{fixed}{Logical vector of length q-1.}
}
\details{
A full q by q Toeplitz correlation matrix is parameterized by q-1 reflection coefficients.  Each PACF is mapped from an unconstrained working coordinate with \code{tanh}; the reflection coefficients are converted to a stable autoregressive representation and then to the implied Toeplitz correlation sequence.  This parameterization preserves positive definiteness while allowing every lag correlation to vary indirectly.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(toeplitzm(time), ism(id))
}
}
\seealso{
\code{\link{ar1m}}, \code{\link{ar2m}}, \code{\link{ar3m}}, \code{\link{mam}}, \code{\link{vsm}}.
}

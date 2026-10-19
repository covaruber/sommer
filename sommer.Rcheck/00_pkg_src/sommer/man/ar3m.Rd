\name{ar3m}
\alias{ar3m}
\title{Third-order autoregressive covariance structure}
\description{
Third-order autoregressive covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
ar3m(x, pacf = c(0.20, 0.10, 0.05), fixed = NULL,
     variance = c("homogeneous", "heterogeneous"), values = NULL)
}
\arguments{
  \item{x}{Ordered factor defining covariance levels.}
  \item{pacf}{Length-three vector of starting partial autocorrelations, each strictly between -1 and 1.}
  \item{fixed}{For homogeneous variance, a logical vector of length three. For heterogeneous variance, a logical vector of length q+2: three PACFs followed by q-1 variance-ratio parameters.}
  \item{variance}{Whether the AR correlation has homogeneous or heterogeneous marginal variances.}
  \item{values}{Optional q positive starting variances when \code{variance = "heterogeneous"}.}
}
\details{
A stationary AR(3) covariance is parameterized through three partial autocorrelations. The PACF parameterization guarantees stationarity while allowing unconstrained working coordinates. Heterogeneous mode uses \deqn{K=D R D} and estimates q-1 positive variance ratios in addition to the PACFs.
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

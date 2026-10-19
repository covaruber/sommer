\name{mam}
\alias{mam}
\alias{ma1m}
\alias{ma2m}
\title{Moving-average covariance structures of order one or two}
\description{
Moving-average covariance structures of order one or two. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
mam(x, order = 1L, theta = NULL, fixed = NULL)
ma1m(x, theta = 0.15, fixed = FALSE)
ma2m(x, theta = c(0.15, 0.05), fixed = NULL)
}
\arguments{
  \item{x}{Ordered factor defining covariance levels.}
  \item{order}{Moving-average order; currently 1 or 2.}
  \item{theta}{Starting MA polynomial coefficients.  Defaults to 0.15 for each order.}
  \item{fixed}{Logical vector with one value per MA coefficient.}
}
\details{
Using the conventional polynomial \eqn{e_t+\theta_1e_{t-1}+\theta_2e_{t-2}}, the autocovariance at lag \eqn{h} is proportional to \eqn{\sum_{j=0}^{p-h}\theta_j\theta_{j+h}} with \eqn{\theta_0=1}.  The covariance is normalized to correlation scale and is zero beyond lag \eqn{p}.  Invertibility of the MA polynomial is not required merely to define this covariance matrix.  \code{ma1m} and \code{ma2m} are convenience wrappers.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(ma1m(time), ism(id))
vsm(ma2m(time), ism(id))
}
}
\seealso{
\code{\link{ar1m}}, \code{\link{toeplitzm}}, \code{\link{vsm}}.
}

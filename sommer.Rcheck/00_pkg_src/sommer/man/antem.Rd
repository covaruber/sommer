\name{antem}
\alias{antem}
\title{Antedependence covariance structure}
\description{
Antedependence covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
antem(x, order = 1L, beta = NULL, innovations = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Ordered variable defining q covariance levels.}
  \item{order}{Antedependence order, 1 <= order < q.}
  \item{beta}{Optional starting regression coefficients for allowed subdiagonals.}
  \item{innovations}{Optional q positive innovation variances.}
  \item{fixed}{Logical vector controlling all beta coefficients and q-1 innovation-variance ratios.}
}
\details{
The modified-Cholesky representation is \deqn{Ty=e,\qquad \mathrm{Cov}(e)=D,\qquad K=T^{-1}DT^{-\mathsf T}.} The matrix \eqn{T} is unit lower triangular, with nonzero regression coefficients only within the requested number of preceding levels.  The first innovation variance is the scale reference; the remaining innovation variances are positive ratios.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(antem(time, order=2), ism(id))
}
}
\seealso{
\code{\link{ar1m}}, \code{\link{toeplitzm}}, \code{\link{vsm}}.
}

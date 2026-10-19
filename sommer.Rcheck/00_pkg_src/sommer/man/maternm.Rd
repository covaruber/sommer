\name{maternm}
\alias{maternm}
\title{Matern spatial covariance structure}
\description{
Matern spatial covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
maternm(x, range = NULL, nu = 0.5, fixed = c(FALSE, FALSE),
        distance = NULL, anisotropy = c("none", "geometric"),
        angle = pi/2, ratio = 1.5)
}
\arguments{
  \item{x}{Numeric coordinate vector or matrix/data frame whose rows give spatial coordinates.}
  \item{range}{Positive starting range parameter.  If NULL, the median positive pairwise distance is used.}
  \item{nu}{Positive starting Matern smoothness parameter.}
  \item{fixed}{Logical of length one or two (range, nu), or four with geometric anisotropy (range, nu, angle, ratio); length one is recycled.}
  \item{distance}{Optional finite symmetric q by q distance matrix with zero diagonal, aligned to the unique spatial locations.}
  \item{anisotropy}{\code{"geometric"} rotates two-dimensional coordinates by \code{angle} and divides the second rotated axis by \code{ratio} before computing distances (ASReml anisotropic Matern).}
  \item{angle}{Starting anisotropy angle in (0, pi).}
  \item{ratio}{Starting positive anisotropy ratio.}
}
\details{
For distance \eqn{d}, the correlation is \deqn{K(d)=\frac{2^{1-\nu}}{\Gamma(\nu)}z^{\nu}K_{\nu}(z),\qquad z=\frac{\sqrt{2\nu}d}{\phi},} where \eqn{\phi} is \code{range}.  Both range and smoothness are positive and are optimized on log scales.  The implementation uses the generic R evaluator and central factor-level numerical derivatives.  Repeated coordinate rows are mapped to the same covariance level.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(maternm(cbind(xcoord,ycoord)), ism(genotype))
vsm(maternm(time), ism(id))
}
}
\seealso{
\code{\link{metricm}}, \code{\link{sar}}, \code{\link{car}}, \code{\link{vsm}}.
}

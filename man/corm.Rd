\name{corm}
\alias{corm}
\title{Uniform correlation covariance shape}
\description{
Uniform correlation covariance shape. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
corm(x, rho = 0.10, fixed = FALSE)
}
\arguments{
  \item{x}{Variable defining covariance levels.}
  \item{rho}{Starting common correlation.}
  \item{fixed}{Logical indicating whether rho is fixed.}
}
\details{
Under the one-scale \code{vsm} architecture, \code{corm} has the same covariance shape as \code{csm}: a unit diagonal and one common off-diagonal correlation.  The overall homogeneous variance is supplied by \code{vsm}.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(corm(environment), ism(genotype))
}
}
\seealso{
\code{\link{csm}}, \code{\link{corvm}}, \code{\link{corhm}}, \code{\link{vsm}}.
}

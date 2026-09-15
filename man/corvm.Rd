\name{corvm}
\alias{corvm}
\title{Uniform correlation covariance shape with corv model label}
\description{
Uniform correlation covariance shape with corv model label. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
corvm(x, rho = 0.10, fixed = FALSE)
}
\arguments{
  \item{x}{Variable defining covariance levels.}
  \item{rho}{Starting common correlation.}
  \item{fixed}{Logical indicating whether rho is fixed.}
}
\details{
\code{corvm} returns the same estimable covariance shape as \code{csm}/\code{corm}.  In older parameterizations a CORV structure may include an additional homogeneous variance, but the current architecture assigns exactly one overall variance to \code{vsm}; adding another scale would be unidentifiable.  The model label is retained for interpretation.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(corvm(environment), ism(genotype))
}
}
\seealso{
\code{\link{csm}}, \code{\link{corm}}, \code{\link{corhm}}, \code{\link{vsm}}.
}

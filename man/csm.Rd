\name{csm}
\alias{csm}
\title{Compound-symmetry covariance structure}
\description{
Compound-symmetry covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
csm(x, rho = 0.10, fixed = FALSE)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{rho}{Starting common off-diagonal correlation.}
  \item{fixed}{Logical indicating whether rho is fixed.}
}
\details{
\code{csm} uses unit diagonal and a common off-diagonal correlation, \deqn{K_{ii}=1, \quad K_{ij}=\rho\ (i\ne j).} Positive definiteness requires \eqn{-1/(q-1)<\rho<1}.  The optimizer uses a bounded-logit working coordinate that maps exactly to this interval.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(csm(environment), ism(genotype))
}
}
\seealso{
\code{\link{corm}}, \code{\link{corvm}}, \code{\link{corhm}}, \code{\link{vsm}}.
}

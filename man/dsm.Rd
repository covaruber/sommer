\name{dsm}
\alias{dsm}
\title{Diagonal heterogeneous covariance structure}
\description{
Diagonal heterogeneous covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
dsm(x, values = NULL, fixed = NULL, theta = NULL)
}
\arguments{
  \item{x}{Variable or design defining the covariance levels.}
  \item{values}{Optional positive starting variances, one per level.}
  \item{fixed}{Logical vector of length q-1 indicating which variance-ratio parameters are fixed.}
  \item{theta}{Optional q by q diagonal starting covariance matrix; if supplied, its diagonal replaces values.}
}
\details{
For \eqn{q} levels, \code{dsm} uses the first variance as the internal scale reference and estimates \eqn{q-1} positive variance ratios, \deqn{K=\mathrm{diag}(1,r_2,\ldots,r_q).} The ratios are optimized on log scales and reported on their natural positive scale.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(dsm(environment), ism(genotype))
}
}
\seealso{
\code{\link{vsm}}, \code{\link{mmes}}.
}

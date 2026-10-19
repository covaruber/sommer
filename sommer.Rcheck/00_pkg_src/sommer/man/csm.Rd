\name{csm}
\alias{csm}
\title{Compound-symmetry covariance structure}
\description{
Compound-symmetry covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
csm(x, rho = 0.10, fixed = FALSE,
    variance = c("homogeneous", "heterogeneous"), values = NULL)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{rho}{Starting common off-diagonal correlation.}
  \item{fixed}{For homogeneous variance, a logical indicating whether rho is fixed. For heterogeneous variance, a logical vector of length q: rho followed by q-1 variance-ratio parameters.}
  \item{variance}{Whether the correlation has homogeneous or heterogeneous marginal variances.}
  \item{values}{Optional q positive starting variances when \code{variance = "heterogeneous"}.}
}
\details{
With \code{variance = "homogeneous"}, \code{csm} uses unit diagonal and a common off-diagonal correlation, \deqn{K_{ii}=1, \quad K_{ij}=\rho\ (i\ne j).} With \code{variance = "heterogeneous"}, it uses \deqn{K=D C(\rho)D,} where the first variance is the reference and the remaining q-1 variances are estimated as positive ratios. Positive definiteness requires \eqn{-1/(q-1)<\rho<1}. The product-level variance remains owned by \code{vsm}.
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
\code{\link{corgm}}, \code{\link{dsm}}, \code{\link{vsm}}.
}

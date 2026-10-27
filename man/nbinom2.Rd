\name{nbinom2}
\alias{nbinom2}
\title{Negative-binomial NB2 family with a log link}
\description{
Creates a sommer-native negative-binomial family with variance
\eqn{\mu+\mu^2/\theta}. It can be used in native Laplace fits without another
package providing the family or fitting engine.
}
\usage{
nbinom2(theta = 1, link = "log")
}
\arguments{
  \item{theta}{Positive starting shape. Larger values approach a Poisson distribution.}
  \item{link}{Currently \code{"log"} only.}
}
\details{
For \code{mmes()} with \code{task$glmm$method="laplace"}, shape is estimated
jointly with fixed effects and free random covariance coordinates by default.
Use \code{task$glmm$familyParameters=list(theta=2, estimate=FALSE)} to fix shape.
For \code{familym()}, nest this configuration under the trait's name.
The constructor's \code{theta} is the default starting value; an explicitly
supplied control value overrides it. Counts must be non-negative integers.

The full likelihood includes the negative-binomial normalization constants.
Native derivatives include shape dependence of both the conditional mode and
Laplace determinant correction. The fitted shape and approximate uncertainty
are available as \code{theta.nb}, \code{theta.nb.se}, and
\code{familyParameters} for a single-family fit. Uncertainty requires a
converged, interior fit with valid joint curvature. Shape near the Poisson
limit can be weakly identified; numerical boundary/curvature warnings are not
calibrated hypothesis tests.

This family also supplies the ordinary R family methods needed by PQL, where
shape is fixed unless the separate PQL shape-estimation control is enabled.
PQL is not changed into a marginal-likelihood estimator by choosing this family.
}
\value{An object of class \code{"family"}.}
\seealso{\code{\link{mmes}}, \code{\link{familym}}, \code{\link{mmes_control}}}
\examples{
family <- nbinom2(theta=2)
family$variance(c(1, 2, 3))
control <- mmes_control(task=list(objective="ml",
  glmm=list(method="laplace", familyParameters=list(theta=2, estimate=FALSE))))
}
\keyword{models}

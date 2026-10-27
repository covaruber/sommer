\name{zipoisson}
\alias{zipoisson}
\title{Zero-inflated Poisson family for native Laplace models}
\description{
Constructs a zero-inflated Poisson family with a Poisson log-mean component
and a separate fixed-effects logit component for structural zeros. An
intercept-only probability is estimated by default by native Laplace fitting.
}
\usage{
zipoisson(pi = 0.1, link = "log", zero = ~1)
}
\arguments{
  \item{pi}{Starting structural-zero probability, a finite value strictly
  between zero and one.}
  \item{link}{Link for the Poisson count mean; currently only \code{"log"}
  is supported.}
  \item{zero}{One-sided fixed-effects formula for the logit structural-zero
  probability. Defaults to \code{~1}, retaining the scalar intercept-only
  interface. Covariates, factors, transformations and offsets are supported;
  random effects are not. A nontrivial formula requires \code{data} and
  \code{estimate=TRUE}.}
}
\details{
For observation \eqn{i}, let \eqn{\mu_i = \exp(\eta_i)} be the Poisson count
mean and \eqn{\pi_i} the structural-zero probability. The response is zero with
probability \eqn{\pi_i + (1-\pi_i)\exp(-\mu_i)}, and a positive count \eqn{y_i} has
probability \eqn{(1-\pi_i)\operatorname{Poisson}(y_i;\mu_i)}.

For \code{zero=~1}, use
\code{familyParameters = list(zeroInflation = value, estimate = FALSE)} in
\code{task$glmm} to fix its probability. For a nontrivial formula,
\eqn{\operatorname{logit}(\pi_i)=w_i'\gamma+a_i}, with fixed design row
\eqn{w_i} and optional zero-component offset \eqn{a_i}. The intercept starts
at \code{qlogis(pi)} and other coefficients at zero; without an intercept,
all coefficients start at zero. Formula coefficients are jointly estimated
and reported in \code{fit$zeroCoefficients} on the logit scale, not as one
scalar probability in \code{familyParameters}. \code{summary(fit)} adds
post-fit joint-curvature standard errors. Zero-only missing covariates use
\code{naMethodX}; in \code{familym()} they affect only their ZIP trait's rows.
The family is currently supported
only by native Laplace ML fitting, including as a named component of
\code{\link{familym}}. PQL fitting is not supported.

The link-scale predictor describes the underlying Poisson count mean. Fitted
and response-scale predicted values are the mixture mean
\eqn{(1-\pi_i)\mu_i}. Simulation draws the structural-zero indicator and then,
when not structurally zero, a Poisson count. Randomized quantile residuals
use the zero-inflated response CDF.

Random effects in the zero component, hurdle models and zero-inflated NB2
are not supported. Sparse data and probabilities near zero or one can yield
weak identification, unavailable uncertainty or numerical-boundary errors.
Convergence alone does not establish a structural-zero mechanism.
The fixed-effects formula is a bounded extension: development calibration
supported moderate-rate interior settings but did not meet its registered
convergence, coefficient-bias and uncertainty-availability criteria in a
low-rate setting. It is not a universally validated mixture-inference method.
}
\value{An object inheriting from \code{"family"}.}
\seealso{\code{\link{mmes}}, \code{\link{familym}}, \code{\link{nbinom2}}}
\examples{
zipoisson(pi = 0.2)
zipoisson(pi = 0.2, zero = ~ treatment)
}
\keyword{models}

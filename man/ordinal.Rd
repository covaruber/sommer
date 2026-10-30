\name{ordinal}
\alias{ordinal}
\title{Ordered categorical family specification}
\description{
Creates a cumulative-logit or cumulative-probit ordered-response family for
single-family native Laplace models fitted with \code{\link{mmes}}.
This experimental implementation requires maximum likelihood and the
Henderson formulation; it is not supported by PQL or \code{familym()}.
}
\usage{
ordinal(link = "logit", thresholds = NULL)
}
\arguments{
  \item{link}{Character string naming the cumulative link. Supported values are
  \code{"logit"} and \code{"probit"}.}
  \item{thresholds}{Optional strictly increasing numeric threshold values used
  to initialize the ordered category cutpoints. To hold them fixed, set
  \code{task$glmm$familyParameters$estimate=FALSE}.}
}
\details{
For categories \eqn{1,\ldots,K}, the cumulative-link convention is
\eqn{P(Y_i \leq k \mid u)=F(c_k-\eta_i)}, with
\eqn{\eta_i=x_i'\beta+z_i'u+o_i}. Here \eqn{F} is the standard logistic or
standard normal CDF and \eqn{c_1<\cdots<c_{K-1}} are cutpoints.
Thus positive predictor coefficients shift probability toward higher
categories. The internal kernel evaluates log category probabilities,
expected category codes, and predictor derivatives through fourth order for
supplied fixed cutpoints. Expected codes are not category probabilities and
depend on the numerical category coding.

The native objective jointly estimates a free first cutpoint and
log-positive successive gaps along with fixed effects and random covariance
parameters. Analytic cutpoint scores include random-mode and Laplace
determinant dependence. It rejects fixed-design confounding of cutpoint
locations and requires observations in every declared category when estimating
cutpoints. Response preparation preserves explicitly ordered factor levels,
or validates numeric integer category codes \eqn{1,\ldots,K} with no gaps.
The response schema is checked before row filtering: unordered factors,
characters, matrices and non-integer codes are rejected rather than sorted,
collapsed or rounded. Ordered factors can declare absent levels, but estimating
their cutpoints requires observations in every declared category after the
standard response, fixed-covariate and random-covariate missing-data policies
are applied. With fixed cutpoints, absent ordered-factor levels are preserved.

Preparation reads \code{thresholds} and \code{estimate} from the
existing \code{task$glmm$familyParameters} list. Control-supplied thresholds
override constructor thresholds. With \code{estimate=TRUE} (the default),
unsupplied cutpoints start at the link-transformed retained cumulative category
frequencies. With \code{estimate=FALSE}, cutpoints must be supplied.
In both cases exactly \eqn{K-1} cutpoints are required.
Estimating cutpoints removes the fixed-design intercept column without
rebuilding the factor contrasts; equivalent constant columns still error.
Fixed cutpoints retain the formula's intercept. Thus a no-intercept formula
that produces a full set of factor indicators is still confounded with free
cutpoint locations and is rejected.

An unrestricted predictor intercept and unrestricted cutpoint locations cannot
both be identified.

The fitted object stores \code{thresholds} (cutpoint estimates with initially
missing standard errors), \code{ordinal$levels} and
\code{category.probabilities}. By default \code{fitted(fit)} and
\code{predict(fit)} return an observation-by-category probability matrix
with category labels as column names. Types \code{"response"} and
\code{"probabilities"} are equivalent; \code{"link"} returns the latent
predictor and \code{"mean"} returns the expected numeric category code.
The internal \code{fitted.values} vector also contains expected codes, not
category probabilities. Such means and response, Pearson and working
residuals depend on the chosen category coding. Deviance residuals use the
observed category log probability; randomized quantile residuals use its
cumulative-probability interval.

Conditional prediction uses fitted modes for known random levels.
\code{re.form=NA} sets random effects to zero. Population prediction
(\code{re.form="population"}) integrates category probabilities over the
supported random covariance structures, using exact normal convolution for
probit and 24-node Gaussian quadrature for logit. Public simulations return
ordered factors with the original category levels; bootstrap refits use
numeric category codes internally and report absent-category failures.

All uncertainty is post-fit: \code{summary()} reports cutpoint standard
errors after transforming the joint marginal Hessian from log-gap coordinates,
and \code{vcov()} retains the fixed-effect covariance convention.
\code{predict(se.fit=TRUE)} includes joint cutpoint, fixed-effect and covariance
parameter uncertainty; conditional prediction also includes random-mode
posterior uncertainty. Probability standard errors have the same shape as
the probability matrix. With \code{cov.fit=TRUE}, the covariance is ordered
column-major: all rows for the first category, then all rows for the second,
and so on. \code{postPEV()} propagates joint parameter uncertainty through
the modes rather than silently conditioning on estimated cutpoints.
Nonconverged, near-boundary or non-positive-definite fits do not receive
success-shaped standard errors. No fit-time numerical Hessian or PEV is
computed, and no new package dependency is required.

Implementation tests do not establish coverage calibration or qualification
for sparse categories, separation or weakly identified covariance parameters.
An external targeted study of 3,000 paired datasets (both links, four
categories and random intercepts) found all 1,200 regular-scenario fits
converged and parameter-recovery gates passed, with 99.7--100 percent joint
SE availability. Nevertheless 18 of 712 registered gates failed, mainly
pointwise probability coverage/RMSE, so the full targeted inference envelope
is not qualified. Regular probability coverage ranged from 90.3 to 98.0
percent; some Wilson intervals excluded nominal 95 percent.
Sparse scenarios intentionally rejected absent categories and weak-variance
scenarios often lacked usable uncertainty despite convergence.
These findings do not qualify random-effect PEV coverage, fixed-cutpoint
inference or general covariance. These findings motivated a further
post-fit probability-interval calibration study.
A subsequent held-out study of 1,200 fresh datasets compared paired Wald,
logit-, probit- and arcsine-transformed probability intervals with fixed
acceptance gates. They failed 8, 3, 4 and 5 of 96 pointwise cell gates
respectively. None qualified across the full envelope, despite complete fit
convergence and 99.7--100 percent SE availability. Transformations alone
are not a demonstrated general coverage correction; the public SE API
remains unchanged. These were delta intervals, not a refit-bootstrap study.
The scalar link functions in this specification describe the cumulative link,
not the expected category code. Scalar variance, deviance, AIC and GLM
initialization methods explicitly error rather than returning placeholder
results. The constructor itself does not validate any response.
}
\value{An object of class \code{"family"}.}
\seealso{\code{\link{mmes}}, \code{\link{familym}}, \code{\link{zipoisson}}}
\examples{
ordinal()
ordinal(link = "probit", thresholds = c(-1, 0.5, 2))
\dontrun{
ctrl <- mmes_control(task=list(objective="ml", glmm=list(method="laplace")))
fit <- mmes(rating ~ treatment, random=~group, data=dat,
  family=ordinal(), control=ctrl)
fitted(fit)                       # labelled category probabilities
predict(fit, re.form=NA)
summary(fit)$thresholds           # uncertainty computed after fitting
simulate(fit, nsim=2, seed=10)
}
}
\keyword{models}

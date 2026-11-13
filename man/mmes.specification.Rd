\name{mmes.specification}
\alias{mmes.specification}
\alias{coef.mmes.specification}
\alias{fitted.mmes.specification}
\alias{residuals.mmes.specification}
\alias{predict.mmes.specification}
\alias{simulate.mmes.specification}
\alias{summary.mmes.specification}
\alias{vcov.mmes.specification}
\alias{confint.mmes.specification}
\alias{anova.mmes.specification}
\alias{vcDelta.mmes.specification}
\alias{logLik.mmes.specification}
\alias{nobs.mmes.specification}
\alias{AIC.mmes.specification}
\alias{BIC.mmes.specification}
\alias{print.mmes.specification}
\alias{print.summary.mmes.specification}
\title{Results for bounded formula specification models}
\description{
Point outputs, simulation and explicitly requested local inference for
formula expressions fitted by the bounded dense reference engine.
}
\usage{
\method{coef}{mmes.specification}(object, ...,
  type=c("mean", "all", "working", "group"))
\method{fitted}{mmes.specification}(object, type=c("response", "link"), ...)
\method{residuals}{mmes.specification}(object,
  type=c("response", "pearson"), ...)
\method{predict}{mmes.specification}(object, newdata=NULL,
  type=c("response", "link"), target=c("conditional", "zero", "population"),
  se.fit=FALSE, cov.fit=FALSE,
  interval=c("none", "confidence", "prediction"), level=.95, trials=NULL, ...)
\method{simulate}{mmes.specification}(object, nsim=1, seed=NULL,
  target=c("conditional", "zero", "population"), ...)
\method{summary}{mmes.specification}(object, ..., inference=FALSE)
\method{vcov}{mmes.specification}(object, ...,
  type=c("mean", "all", "working", "random", "conditional", "group"))
\method{confint}{mmes.specification}(object, parm, level=.95, ...,
  type=c("mean", "all", "working"))
}
\arguments{
\item{object}{A fitted \code{mmes.specification} model.}
\item{type}{Output scale or parameter/covariance layout.}
\item{newdata}{Optional nonempty prediction data frame. Random levels must be
known for conditional and population targets.}
\item{target}{Conditional fitted modes, zero deviations, or the integrated
prior population. Population link prediction integrates the nonlinear link
predictor, not its zero-deviation value.}
\item{se.fit, cov.fit}{Request standard errors and, optionally, full prediction
covariance ordered as \code{unlist(fit, use.names=FALSE)}.}
\item{interval}{No interval, local normal confidence interval, or local normal
future-response interval. The latter requires scalar response output.}
\item{level}{Confidence level strictly between zero and one.}
\item{trials}{Positive integer trials for binomial prediction intervals; one
or one per prediction row. Required for new-data grouped-binomial intervals.}
\item{nsim}{Positive integer simulation count.}
\item{seed}{Optional seed. Seeded calls restore the caller's RNG state.}
\item{inference}{Explicitly request coefficient uncertainty in the summary.}
\item{parm}{Unique parameter names or indices. Omitted selects all.}
\item{...}{No additional arguments are accepted.}
}
\details{
Default coefficients are fixed regression coefficients on their declared
working/link scales. \code{type="working"} includes every free optimizer
coordinate; \code{"all"} reports mean coefficients and natural covariance and
family nuisance coordinates, including fixed entries with zero variance.
\code{coef(type="group")} returns natural parameter values per retained
predictor row, including fitted parameter deviations.

Scalar response predictions are vectors. A single ordinal or nominal family
returns a probability matrix with declared category ordering; mixed
categorical families return trait-named lists. Link output preserves nominal
logit matrices. Categorical response residuals subtract fitted category
indicators (nominal) or expected category codes (ordinal).

Population integration uses the joint covariance of all active predictor
coordinates, including correlated nominal logits within a record. Separate
parameter records remain independent. Quadrature doubles order until the
requested tolerance or fails explicitly at its order/memory limit.
Simulation draws the fitted conditional modes, zero deviations or independent
declared prior covariance blocks before applying parameter links. Binomial
predictions and simulated binomial responses are proportions, using the
declared integer trials.

Conditional curvature composes symbolic expression derivatives with the full
native family curvature. Outer scores use
central numerical differences. Legacy noncanonical scalar links use a
numerical second derivative of the inverse link.

\code{post_nonlinear()} and covariance/interval requests verify reconstruction,
convergence, stationarity, numerical derivative agreement, boundaries and
positive identified curvature. Conditional modes are re-certified at tighter
tolerances for post-fit derivatives, without changing the saved fitting controls
or statistical model. Singular information is not repaired using a
pseudoinverse. \code{vcov(type="conditional")} holds all marginal
parameters fixed; \code{"random"} adds mode sensitivity and is a local
empirical Bayes PEV approximation, not an exact frequentist PEV.
\code{"random"}, \code{"conditional"}, and \code{"group"} are unavailable
for ordinary fixed-only fits.
\code{"group"} propagates parameter links and mode motion to the natural
group parameters. No uncertainty is computed during fitting or ordinary
summary. All interval results are local normal approximations, not exact
discrete or simultaneous intervals. Probability confidence limits may extend
outside zero and one. Population prediction intervals include both integrated
family variance and prior variation.

Ordinary fixed-only models retain their \code{lm}/\code{glm} covariance
semantics, with fitted ML dispersion when appropriate. For nonlinear models,
all estimated covariance and nuisance coordinates enter the checked observed
information. \code{vcDelta()} names the non-mean natural reporting coordinates
\code{V1}, \code{V2}, and so on, and returns their mapping as an attribute.
Comparisons admit only regular nested ordinary fixed terms with unchanged
nonlinear, family, random and observation semantics; Wald tests require named
mean coordinates and explicit finite nulls. Boundary-variance tests are not
calibrated.
}
\value{
Coefficient vectors or matrices, predictions in the fitted family layout,
simulation data frames (trait lists for mixed categorical responses),
covariance matrices, confidence intervals, or a summary. Prediction uncertainty
returns \code{fit}, \code{se.fit}, an approximation label, optional
\code{cov.fit}, and optional \code{lower}/\code{upper} with matching shapes.
}
\seealso{\code{\link{parameter}}, \code{\link{nonlinear_curve}},
\code{\link{post_nonlinear}}, \code{\link{profile_nonlinear}}}

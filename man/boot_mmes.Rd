\name{boot_mmes}
\alias{boot_mmes}
\title{Post-fit parametric bootstrap for native Laplace models}
\description{
Simulates new population random effects and responses and refits the native
Laplace likelihood without calculating PEV or numerical-Hessian uncertainty.
Supports parameter-recovery summaries and simulation-based comparisons of nested
fixed-effect models.
}
\usage{
boot_mmes(object, object2 = NULL, nsim = 199L, seed = NULL, maxit = 1000L)
}
\arguments{
\item{object}{A converged native Laplace \code{mmes} fit. For a comparison,
the smaller, generating null model.}
\item{object2}{Optional converged larger model differing only in strictly
nested count fixed effects. Response rows, offsets, families, trials,
zero designs and random covariance models must match.}
\item{nsim}{Positive integer number of planned bootstrap datasets.}
\item{seed}{Optional nonnegative integer seed. When supplied, the caller's
random-number state is restored on exit.}
\item{maxit}{Positive integer outer iteration limit for each refit.}
}
\details{
This is an opt-in post-fit calculation, not part of fitting. It uses the
stored retained designs, offsets, trials and native likelihood specifications,
so no evaluation of the original data expression or response formula is needed.
Grouped binomial simulations are converted from proportions back to counts.
Each refit starts at that model's fitted parameters and uses native first
derivatives and the stored mode controls. Fits are sequential; no dependencies
or parallel workers are added. Output storage is subject to the generating
fit's \code{execution$denseMemoryMB} budget.

Single-model percentile intervals are descriptive bootstrap intervals on the
stored optimizer scale: count coefficients, covariance descriptor coordinates,
and log shape/dispersion or zero-logit coefficients. They are not automatically
transformed to natural covariance parameters and are not studentized or
bias-corrected. They use successful refits only and require separate coverage
calibration, especially at boundaries.

Single-family \code{\link{ordinal}} models retain their declared category
labels and cutpoint policy for comparisons. Cutpoint draws use the first
location followed by log-positive successive gaps, not natural cutpoint
coordinates. Simulations use numeric codes internally; refits with an absent
declared category fail explicitly when cutpoints are estimated. Ordinal
comparisons may differ only by nested main-predictor fixed effects, not by
category labels or fixed thresholds.

Single-family \code{\link{multinomial}} fits use their joint softmax
likelihood in bootstrap refits. Simulated factor labels are mapped back to
the fitted internal reference coding; all declared categories remain in the
response schema even when a bootstrap sample omits one.

With two models, simulations come from the smaller model, drawing new random
effects rather than holding conditional modes fixed. The statistic is twice
the difference in approximate marginal log likelihood. If every replicate
succeeds, the p-value is \eqn{(1+e)/(1+B)}, where \eqn{e} is the number of
bootstrap statistics at least as large as observed and \eqn{B} is
\code{nsim}. Monte Carlo uncertainty is reported separately. Use substantially
more than 199 replicates for precise tail probabilities.

Failures and nonconvergence are retained, not retried until success or
silently treated as non-rejections. If any replicate fails, \code{p.value}
is \code{NA}; \code{p.value.bounds} assigns all failures first as
non-exceedances and then as exceedances. These are failure-sensitivity
bounds, not confidence limits. Boundary replicates are retained and flagged.
A boundary generating fit produces a warning; this is not a general
variance-component boundary test, zero-inflation selection test, or a cure
for an unidentified generating model.

Parameter simulations condition on the fitted generating parameters.
No numerical Hessian, PEV, or full fitted objects are produced for refits.
The original fitted objects are not modified. Ordinary Gaussian LMM and
PQL fitting paths are unchanged.
}
\value{
A list with \code{estimates}, optimizer-scale \code{draws},
\code{percentile} intervals, \code{successful} and \code{nsim} counts,
and one row per planned replicate in \code{diagnostics}. Diagnostics include
convergence, likelihoods, score magnitude, rejected evaluations, boundaries,
warnings and errors. Comparisons also return \code{statistic},
\code{p.value}, \code{p.value.bounds} and \code{monteCarloSE}.
}
\seealso{\code{\link{mmes}}, \code{\link{simulate.mmes.glmm.laplace}}}
\examples{
\dontrun{
boot <- boot_mmes(fit, nsim=999, seed=17)
boot$diagnostics
boot$percentile
test <- boot_mmes(smaller, larger, nsim=999, seed=17)
test$p.value
test$p.value.bounds
}
}

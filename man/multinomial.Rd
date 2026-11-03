\name{multinomial}
\alias{multinomial}
\title{Nominal multinomial family specification}
\description{
Creates a nominal, unordered multinomial response family for native Laplace
mixed-model fitting with \code{\link{mmes}}.
}
\usage{
multinomial(reference = NULL)
}
\arguments{
  \item{reference}{Optional label of the category used as the zero-logit
  reference. By default, the last declared response-factor level is used.}
}
\details{
This family supports unordered factors and plain integer category codes
\code{1, ..., K}. Ordered factors, character responses, matrices, noninteger
codes, and numeric codes with gaps are rejected. The response must contain
observations in every declared category after missing-data filtering.

For \eqn{K} categories the model estimates \eqn{K-1} logits against the
reference category. For example, with \code{reference="green"}, the logits
are \eqn{\log(p_{red}/p_{green})} and
\eqn{\log(p_{blue}/p_{green})}. A softmax calculation converts these
comparisons to probabilities that sum to one. Changing the reference changes
the coefficient representation, not the category labels returned by public
methods.

Fit with maximum likelihood and the native Laplace method, for example:
\code{mmes(y ~ x, random = ~ group, data = dat, family = multinomial(),
control = mmes_control(task = list(objective = "ml",
glmm = list(method = "laplace"))))}. At least one random-effect term is
required. The model uses category-specific fixed coefficients. Random effects
can be shared across comparisons (for example \code{~group}) or category
specific using the internal \code{.nominalCategory} factor in the random
formula. The latter is added to the stacked design by the fitter.

The public fit requires the Henderson formulation and \code{output$inverse=0};
uncertainty is deferred until after fitting. PQL, supplied
\code{W}, weights formulas, \code{vcc}, structured residual covariance,
and non-default Gaussian solver/acceleration controls are unsupported.
Random-effect relationships and covariance structures use the existing
\code{mmes_control()} and covariance-descriptor machinery. No new package
dependency is required.

Categorical \code{\link{familym}} supports nominal traits together with
ordinal and scalar traits, using one numeric long-format response column.
That route uses \code{.familyCoordinate}, not \code{.nominalCategory}, for
explicit category-specific random covariance terms, and returns trait-named
lists from \code{fitted()} and \code{predict()}. Existing random coefficients
are shared, not silently expanded into independent effects.

\code{fitted()} and \code{predict()} return a matrix of category probabilities
with the original response levels as columns; \code{type="link"} returns the
\eqn{K-1} logits. Conditional prediction uses fitted modes for known random
levels, and \code{re.form=NA} sets random effects to zero.
\code{re.form="population"} integrates category probabilities over shared or
correlated category-specific unit random intercepts with identity relationship
precision. Up to four categories and multiple independent terms are supported;
new independent group levels are allowed if the fixed design remains valid.
Random slopes, varying/nonunit loadings, nonidentity relationships, averaged
logits, and population standard errors/covariances are explicitly rejected.
Deterministic 24- and 32-point tensor normal quadrature on the active covariance
eigenspace must agree within \eqn{10^{-6}} absolute probability; this
consistency check is not a guaranteed integration error bound. Estimated
workspace respects \code{execution$denseMemoryMB}. Population probabilities
are plug-in means under the fitted covariance, not calibrated uncertainty.
\code{simulate()} returns factors preserving
the original category labels and can draw new population random effects with
\code{re.form="population"}.

Post-fit \code{summary()} and \code{vcov()} calculate fixed-effect uncertainty
from the numerical Hessian of the marginal Laplace objective. \code{postPEV()}
calculates random-effect posterior variances after fitting.
\code{predict(se.fit=TRUE)} uses the joint fixed/random conditional curvature;
its probability-scale standard errors use a delta-method transformation.
These default prediction standard errors condition on the fitted covariance
parameters. Opt-in \code{predict(se.fit=TRUE, uncertainty="full")} includes
covariance-parameter uncertainty and movement of conditional random modes
using the full marginal parameter covariance, plus random-mode posterior
covariance for conditional predictions. Zero-random-effect prediction uses
the marginal fixed-coefficient covariance. This post-fit first-order
approximation does not provide calibrated intervals or future-response noise.
Boundary, nonconverged, or non-positive-definite fits do not
receive success-shaped uncertainty. \code{boot_mmes()} supports nominal
parameter-recovery and nested fixed-effect refits by simulating population
random effects and refitting the joint softmax likelihood. Bootstrap draws
are reported on the fixed-coefficient and native covariance-coordinate
scales. Failed or nonconverged replicates remain in diagnostics, and
percentile intervals use successful refits only. This descriptive bootstrap
has not yet been calibrated as a replacement for Wald intervals.

The fitted object's \code{laplaceDiagnostics} separates the optimizer code
and message, \code{scoreMax} (largest absolute optimizer-coordinate score),
covariance boundary flags, and \code{saturatedRows} (rows whose largest fitted
probability is within \eqn{\sqrt{\epsilon}} of one, where \eqn{\epsilon} is
machine precision). Saturation suggests possible separation or weak
identification; it does not prove separation.
\code{covarianceCondition} reports the largest/smallest covariance eigenvalue
ratio per random term. Ratios exceeding \eqn{\epsilon^{-1/3}} trigger an
ill-conditioning advisory, not a new fit rejection or uncertainty gate.
\code{issues} provides actionable messages printed by \code{summary()}, and
\code{uncertaintyStatus} distinguishes deferred calculations from unavailable
uncertainty. These checks do not compute Hessians or PEV during fitting.
Unstructured nominal covariance precision derivatives use triangular-factor
solves to avoid cancellation near singularity.

Likelihood, AIC and BIC are available for the fitted nominal model, and
nested fixed-effect likelihood comparisons are supported when their response,
random structure and covariance specification match.
\code{residuals(type="deviance")} returns unsigned
\eqn{\sqrt{-2\log p_{i,y_i}}} magnitudes, evaluated stably from conditional
logits. Their squares sum to the conditional response deviance, not the
marginal Laplace likelihood deviance. Response/Pearson residuals retain their
labelled observation-by-category matrices; categories within a row are
dependent. Working/quantile residuals remain unsupported because no natural
nominal category ordering is imposed. These post-fit model checks are
descriptive, not calibrated outlier tests; see \code{\link{residuals.mmes}}.
Statistical calibration
for the nominal model remains a separate validation task.
}
\value{An object of class \code{"family"} with class
\code{"sommer_nominal_family"}.}
\seealso{\code{\link{mmes}}, \code{\link{ordinal}}, \code{\link{familym}},
\code{\link{postPEV}}, \code{\link{predict.mmes}}}
\examples{
multinomial()
multinomial(reference = "control")
\dontrun{
ctrl <- mmes_control(task=list(objective="ml", glmm=list(method="laplace")))
fit <- mmes(color ~ treatment, random=~group, data=dat,
  family=multinomial(reference="green"), control=ctrl)
fitted(fit)                       # labelled category probabilities
predict(fit, re.form=NA)          # random effects set to zero
predict(fit, re.form="population") # eligible population category means
summary(fit)$betas                # post-fit fixed-effect uncertainty
simulate(fit, nsim=2, seed=10)
}
}
\keyword{models}

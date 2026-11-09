\name{post_nonlinear}
\alias{post_nonlinear}
\alias{vcov.mmes.nonlinear}
\alias{confint.mmes.nonlinear}
\title{Post-fit inference for traditional nonlinear mixed models}
\description{
Computes explicitly requested local uncertainty after nonlinear fitting.
No curvature inverse, PEV or parameter uncertainty is computed during fitting.
}
\usage{
post_nonlinear(object, mode=1L, control=NULL)
\method{vcov}{mmes.nonlinear}(object, ...,
  type=c("mean", "all", "working", "random", "conditional", "group"))
\method{confint}{mmes.nonlinear}(object, parm, level=0.95, ...,
  type=c("mean", "all", "working"))
}
\arguments{
\item{object}{A converged, locally identifiable nonlinear \code{mmes} fit.}
\item{mode}{\code{0} computes marginal parameter uncertainty only;
\code{1} additionally computes random-effect diagonal approximate PEVs;
\code{2} retains the complete random-effect approximate PEV and
parameter/mode cross-covariance and natural-scale group-parameter covariance matrices.}
\item{control}{Optional \code{mmes_control()} object overriding post-fit
numerical tolerances and memory budget. Defaults to the fitted control.
Other fitting choices are not changed.}
\item{type}{\code{"mean"} matches natural-scale \code{coef(object)};
\code{"all"} includes covariance and family nuisance parameters, including
fixed parameters with zero variance; \code{"working"} returns all free
optimizer-coordinate covariance; \code{"random"} returns the approximate
PEV; \code{"conditional"} returns random-mode covariance conditional on all
estimated marginal parameters; \code{"group"} returns natural-scale group-curve
parameter covariance, ordered by \code{as.numeric(object$groupParameters)}.}
\item{parm}{Unique parameter names or indices for the selected \code{type}.
Omitted means all parameters of that type.}
\item{level}{Confidence level strictly between zero and one.}
\item{...}{No additional arguments are accepted.}
}
\details{
Let \eqn{a} contain all free marginal parameters, including constrained
covariance coordinates and estimated family parameters. Central differences
of the analytic marginal-likelihood score give the observed information
\eqn{I}. Differences at steps \eqn{h} and \eqn{h/2} must agree. The reported
working covariance is \eqn{V=I^{-1}}. Gaussian additive fits use their exact
Gaussian likelihood; other fits use the fitted Laplace likelihood.
Native conditional curvature is exact, including nonlinear second-derivative
terms. Its inverse, restored from the proper latent coordinates, is
\eqn{S}. With \eqn{D=d\hat u/da}, estimated by checked central differences,
the empirical Bayes approximation is
\deqn{\operatorname{Cov}(a,\hat u)=VD',\qquad
 \operatorname{PEV}(u)\ \mathrel{\approx}\ S+DVD'.}
The first matrix is not a sampling covariance of a joint maximum-likelihood
estimate of random effects. The PEV is a local estimated-parameter
approximation, not an exact unconditional frequentist prediction-error
variance. Both quantities are labeled separately. No prior precision inverse
of singular reduced-rank covariance is used.

Natural-scale covariance and group-curve standard errors use a first-order
delta method, including loading-map and conditional-mode motion.
All cross-parameter terms are retained for independent and joint random
curve specifications. \code{vcov(type="group")} explicitly requests full
post-fit group-parameter covariance and requires random curve parameters;
its diagonal agrees with squared \code{groupParameterSE}. It combines
marginal-parameter uncertainty and conditional latent curvature and is a
local empirical Bayes delta approximation, not an exact frequentist PEV.
\code{confint()} uses normal working-coordinate intervals and exponentiates
positive curve coordinates; it does not test positive parameters against zero.
\code{confint(type="working")} returns all free working-coordinate intervals;
\code{type="all"} includes covariance and family nuisance parameters using
natural-scale symmetric delta intervals (which may cross parameter bounds).
Positive curve intervals remain transformed working-normal intervals.
\code{covparams_mmes_se()} applies the same covariance to the existing
descriptor reporting callbacks. Fixed parameter standard errors are zero.
\code{covmatrix_mmes(se=TRUE)} propagates working covariance directly to
the implied covariance and correlation matrices.
\code{vcDelta()} supports scalar formula transforms of the natural covariance,
residual and family nuisance parameters, ordered as in
\code{coef(object, type="all")} after removing fixed/curve entries.
These entries are named \code{V1}, \code{V2}, and so on; the returned
\code{parameters} attribute records this mapping. For an additive Gaussian
model with one random variance these are random variance and residual
variance, respectively. A ratio is a parameter transform, not automatically
a validated reliability measure for a nonlinear response.

Numerical choices use
\code{task$nonlinear$inference=list(step=1e-4, derivativeTol=5e-3,
curvatureTol=1e-8, maxCondition=1e10)} in \code{mmes_control()}.
Steps are relative to \code{max(1,abs(coordinate))}. Positive-definiteness and
conditioning checks use diagonally standardized information.
Nonconverged, nonstationary, weakly identified, indefinite-curvature and
flagged-boundary fits fail explicitly. Structural rank reduction alone is
not treated as an estimated boundary. Workspace is checked against
\code{execution$denseMemoryMB}. The current implementation reconstructs
operation-local native contexts and stores a full conditional covariance;
even mode 1 can require quadratic workspace. No serialized pointers are retained.

\code{predict(..., se.fit=TRUE, cov.fit=TRUE)} returns \code{fit},
\code{se.fit}, \code{cov.fit} and an approximation label. Shapes of fitted
vectors, probability matrices and trait lists are preserved.
\code{cov.fit} orders coordinates as \code{unlist(fit, use.names=FALSE)}
(column-major for matrices). Conditional prediction propagates both
\eqn{S} and total parameter derivatives, including mode motion.
Zero-deviation and prior-population prediction propagate marginal parameter
uncertainty only: prior variability is not an uncertainty of the population
mean. Population derivatives also include estimated covariance and
family-parameter effects through the existing checked quadrature.

\code{predict(..., interval="confidence", level=0.95)} adds \code{lower}
and \code{upper} using a normal first-order approximation on the requested
scale. Probability-scale intervals can exceed zero or one; these are not
simultaneous or exact bounded intervals.
\code{interval="prediction"} additionally includes future-observation
variance, and for population predictions integrates both conditional
variance and squared mean over the joint prior. These are normal
approximations, not exact discrete quantiles; count intervals can extend
below zero. Prediction intervals require scalar response-scale outputs.
Categorical probabilities support confidence intervals, not numeric
future-category intervals. Grouped-binomial new-data prediction intervals
require \code{trials}, one positive integer or one per response. All
binomial outputs are proportions, not success counts.

\code{summary(object)} displays saved uncertainty only;
\code{summary(object, inference=TRUE)} explicitly requests it.
\code{coef(object, type="working")} returns free working coordinates;
\code{type="all"} returns natural-scale parameter values.
\code{AIC()} and \code{BIC()} count all free optimizer coordinates.
\code{anova(smaller, larger)} supports only strictly nested ordinary fixed
terms with identical data, inputs, nonlinear curves, random structure,
constraints and likelihood semantics, at regular identifiable interior fits.
It uses an asymptotic likelihood-ratio reference, approximate for Laplace fits.
\code{anova(object, parm=..., null=...)} requires explicit finite nulls on
working coordinates and returns an approximate Wald test for fixed/curve
parameters only. Automatic nonlinear-term selection, random-variance
boundary tests, \code{wald_mmes()} denominator degrees of freedom,
\code{postVarU()} and linear-model reliability formulas remain unsupported.
No boundary test is silently replaced by an ordinary chi-square test.
}
\value{
\code{post_nonlinear()} returns the fit with a serializable
\code{nonlinearInference} list containing estimates, working/natural
covariance, observed information, numerical settings and method labels.
Random modes additionally populate \code{uPevList}, \code{uPevMode},
conditional covariance, mode sensitivities and group-parameter standard errors.
No linear coefficient inverse \code{Ci} is manufactured.
\code{vcov()} and \code{confint()} return matrices.
Calling these methods does not mutate the input fit.
}
\seealso{\code{\link{nonlinear_curve}}, \code{\link{postPEV}},
\code{\link{mmes_control}}}
\examples{
set.seed(12)
dat <- expand.grid(x=seq(0, 5, length.out=10), group=factor(seq_len(12)))
dat$y <- 3*exp(-0.6*dat$x)+rnorm(12, sd=0.2)[dat$group]+
  rnorm(nrow(dat), sd=0.1)
fit <- mmes(y~0, random=~group, data=dat, nonlinear=nonlinear_curve("x"),
  control=mmes_control(task=list(objective="ml"),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
post <- post_nonlinear(fit)
confint(post)
predict(post, newdata=dat[1:3, ], interval="confidence")
}

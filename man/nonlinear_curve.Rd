\name{nonlinear_curve}
\alias{nonlinear_curve}
\alias{coef.mmes.nonlinear}
\alias{fitted.mmes.nonlinear}
\alias{residuals.mmes.nonlinear}
\alias{predict.mmes.nonlinear}
\alias{simulate.mmes.nonlinear}
\alias{summary.mmes.nonlinear}
\alias{logLik.mmes.nonlinear}
\alias{nobs.mmes.nonlinear}
\alias{AIC.mmes.nonlinear}
\alias{BIC.mmes.nonlinear}
\alias{anova.mmes.nonlinear}
\alias{print.mmes.nonlinear}
\alias{print.summary.mmes.nonlinear}
\title{Traditional nonlinear mixed-model curves}
\description{
Declares a native parametric link-predictor contribution for nonlinear mixed models fitted
by \code{\link{mmes}}.
}
\usage{
nonlinear_curve(input, curve=c("exponential", "logistic", "saturation"),
                random=list(), joint=NULL)
}
\arguments{
\item{input}{Name of one numeric column in \code{data}.}
\item{curve}{One of \code{"exponential"}, \code{"logistic"} or
\code{"saturation"}.}
\item{random}{Empty list or uniquely curve-parameter-named, one-term random
formulas for independent parameter effects, for example
\code{list(rate=~vsm(ism(subject), sigma2=0.04))}. Existing formulas are unchanged.}
\item{joint}{Optional list with \code{parameters} (at least two unique curve
parameter names) and one \code{random} formula, for example
\code{list(parameters=c("asymptote","rate"),
random=~vsm(usm(.parameter),ism(subject)))}. The internal \code{.parameter}
axis has the declared parameter order. Joint and independent parameter
names must be disjoint.}
}
\details{
The supported curves are
\deqn{A\exp(-kx),}
\deqn{A/(1+\exp(-k(x-m))),}
and
\deqn{Ax/(K+x).}
Their named parameters are \code{asymptote}, \code{rate} (exponential and
logistic), \code{midpoint} (logistic), and \code{half} (saturation).
All except \code{midpoint} are strictly positive and optimized on the log scale.
Exponential and saturation inputs must be non-negative.
There is no automatic input scaling.

With \code{random=list(), joint=NULL}, the curve is a fixed mean and ordinary random
effects are additive Gaussian terms. A nonempty \code{random} list declares
Gaussian deviations in named curve coordinates: log scale for positive parameters,
identity for midpoint. For example, random rate means
\eqn{\log(k_j)=\log(k)+b_j}. The shared positive parameter is the
zero-deviation value (and marginal median), not its arithmetic mean.
All curve parameters may vary, optionally together with independent
ordinary additive effects. Separate named formulas specify independent
parameter blocks, even when their grouping variables coincide.
The \code{joint} block permits correlated Gaussian deviations using the
existing \code{vsm()} covariance syntax. It may coexist with independent
formulas for other parameters (for example midpoint).
Cross-covariance with ordinary additive effects remains unsupported.
Each observation must select exactly one unit incidence for each declared
random curve parameter; random slopes and overlapping loadings
are rejected. Shaped native covariance descriptors and known positive-definite
relationship precisions use the existing covariance machinery.
The joint \code{.parameter} axis must occur exactly once as a covariance-shaping
factor, not as the final grouping term. It is generated internally in the
declared order and must not be a data column. No user-data stacking is needed
and no response is replicated. For example \code{usm(.parameter)} specifies
an unstructured covariance among log-asymptote and log-rate deviations;
\code{ism(subject)} supplies the grouping basis and \code{Gu} retains its
usual known-precision meaning. Additional shaping factors and \code{rrm()}
retain their descriptor semantics.
Group-parameter rows align by covariance/group labels across independent
blocks; a parameter with no deviation for a row uses its shared coordinate.
Multiple-parameter covariance term labels are prefixed by
\code{curve[parameter]} (or the comma-separated joint parameter names), so
identical grouping formulas remain distinguishable in reporting and named
constraints. Their parameter indices and existing \code{vcc} mapping semantics
are unchanged. Multiple-parameter \code{curveRandomMode} names include both
the parameter and group-coordinate label.
Joint second and third derivatives include all cross-parameter terms.
Random covariance factors must have native or fixed evaluators and native
analytic (or parameter-free) derivatives. The built-in \code{rrm()} is also
supported through its native normalized loading map, rather than its R
covariance evaluator. Arbitrary R callbacks and numerical covariance
derivatives are rejected. Dedicated numerical checks cover \code{ism()},
\code{dsm()}, \code{usm()}, \code{ar1m()}, \code{csm()}, heterogeneous and
homogeneous \code{fam()}, and \code{rrm()}, including representative Kronecker
products. This is not a blanket guarantee for every covariance constructor.
For Gaussian identity observations, integration is exact Gaussian ML for purely additive
effects and Laplace ML for random curve parameters. Family-aware Laplace
uses the exact conditional Hessian, including likelihood-score-times-second-derivative terms, and full
analytic outer scores including mode motion. This native mixed-model route
requires at least one random term (additive or curve); fixed-only curves use
the reference specification route. Full-coordinate covariance must be positive definite,
except for built-in reduced-rank factors. For \code{rrm()}, the effect is
represented as \eqn{b=L v}, with a proper Gaussian prior on latent \eqn{v}.
Loading parameters remain jointly estimated, and analytic scores include
the moving loading map. No inverse of singular covariance, artificial
specific variance or diagonal ridge is introduced.
Relationship precision and the latent prior must remain positive definite;
exact zero overall variance is not accepted. With additive random effects retained,
remove the curve-random declaration to fit the zero-curve-variance model.
Exact zero is never replaced silently by a
small positive variance. Independent homogeneous Gaussian residuals are supported;
structured observation residuals are not.

Response families reuse the native likelihoods: binomial with logit, probit
or cloglog link; Poisson, NB2, Gamma and zero-inflated Poisson with log link;
Gaussian identity; ordinal logit/probit; nominal softmax; and rowwise
\code{familym()} mixtures, including categorical traits. Apart from the
existing exact Gaussian additive route, these use Laplace ML. The curve
is added to the link predictor, not the response mean. It is shared across
traits and across every nonreference nominal logit; separate nonlinear
parameter sets by category or trait are not introduced. Existing categorical
linear-design conventions still apply, including \code{.nominalCategory}
for a standalone nominal model and \code{.familyCoordinate} for mixtures.
Single categorical responses may be factors; categorical mixtures use
numeric category codes as in ordinary native Laplace models.

Numerical choices belong to \code{\link{mmes_control}}.
\code{optimizer} accepts \code{algorithm="nlminb"} and
\code{gradient="analytic"}. Supply a complete named vector of natural-scale
curve starts through \code{mmes(start=list(parameters=...))}, or omit it
for deterministic automatic initialization.
Supplying \code{start} overrides automatic initialization.
Automatic starts are deterministic and require input variation and a positive
response maximum after subtracting the offset. They do not consume random
numbers. For non-Gaussian models automatic starts use unit link-scale amplitude;
supplied starts are recommended. No silent numerical-gradient or optimizer
fallback is performed. Family nuisance parameters use the existing
\code{familyParameters} list (trait-nested for \code{familym()}),
including fixed/estimated NB2 shape, Gamma/Gaussian dispersion,
zero-inflation probability and ordinal thresholds. Other GLMM fitting controls
do not apply. Zero-component regression uses \code{zipoisson(zero=~...)}.

For Laplace fits, \code{mode} defaults to
\code{list(maxit=100L, tol=1e-8, stepTol=1e-8,
warmStart="auto", starts="deterministic")}. Here \code{"auto"} selects
deterministic cold starts, equivalent to \code{"off"}.
Safeguarded Newton steps compare
zero and \eqn{+/-0.5} curve-coordinate starts (latent curve coordinates for
reduced-rank terms). Final curvature must be
positive definite without damping; failures, score norms, iterations and
competing modes are reported. This is not a global-mode guarantee.
Warm starts and stochastic starts are not supported.

Use \code{convergence$maxIters}, \code{convergence$likelihoodTol} and
\code{convergence$parameterTol} for the outer optimizer. The latter tolerances
map to \code{nlminb} relative objective and parameter-step tolerances.
Resource admission uses \code{execution$denseMemoryMB}; thread policy uses
\code{execution$blasThreads}. Unsupported non-default controls are rejected.
Sparse-support components bound Hessian/factor workspace before sparse
products are assembled; a large connected relationship graph can be refused
conservatively even if its actual fill would fit. This budget is not a
process-wide RSS limit. \code{diagnostics$profile="force"} enables fitting
phase timings, and \code{\link{profile_nonlinear}} measures separate inclusive
kernels, dimensions, work counts and retained object sizes after fitting.

The combined ordinary/nonlinear mean must be locally identifiable.
Starts with rank-deficient mean derivatives are rejected; rank deficiency at
the accepted fit warns. Nonlinear optimization can have local optima;
convergence and final scores are reported, not a guarantee of a global optimum.

Fitting returns point estimates only. PEV and all uncertainty are explicitly
post-fit processes; see \code{\link{post_nonlinear}} for observed parameter
curvature, natural-scale covariance, conditional-mode covariance, approximate
PEV and prediction intervals. \code{postPEV()} dispatches to this nonlinear
implementation; \code{vcov()} and uncertainty prediction arguments request
operation-local post-fit calculations. Summaries display saved uncertainty,
or compute it only when explicitly requested with \code{inference=TRUE}.
\code{predict(re.form=NULL)} uses conditional modes;
\code{predict(re.form=NA)} uses zero random deviations.
\code{predict(re.form="population")} integrates the curve against its fitted
Gaussian prior, not a mode/parameter uncertainty distribution. It generally
differs from zero-deviation prediction for random curves. For non-Gaussian
responses, additive effects also change the population response mean.
Their joint prior, including within-observation nominal predictor covariance,
is integrated; link output integrates the predictor and response output
integrates means or category probabilities. Gaussian quadrature compares increasing orders with
\code{population=list(order=20L, maxOrder=160L, tol=1e-7)};
failure to meet scaled agreement is an explicit error. Conditional and
population new-data predictions require represented random levels and the
fitted covariance basis. Fitted \code{groupParameters} and
\code{curveRandomMode} expose natural-scale group parameters and coordinate modes.
For reduced-rank models, \code{randomMode} and \code{curveRandomMode} are
restored to the original full coordinates; \code{latentRandomMode} and
\code{randomCoordinates} retain the latent point mode and dimension metadata.
Population means use the implied full-coordinate marginal prior variance,
including singular positive-semidefinite reduced-rank covariance.
\code{fitted(type="link")} and \code{predict(type="link")} return link-scale
predictors; response output is numeric for scalar families, category-probability
matrices for ordinal/nominal models, and trait-named lists for categorical
mixtures. Residual types retain the existing native-family conventions;
unordered nominal working/quantile residuals are undefined and rejected.
\code{simulate(nsim=1, seed=NULL, re.form="conditional", maxMemoryMB=NULL)}
draws observations at fitted modes; \code{re.form="population"} draws new
random effects from their fitted prior before evaluating the curve and
response likelihood. Explicit seeds restore the caller's RNG state.
These are generative draws, not PEV or parameter-uncertainty products.
}
\value{
A versioned, serializable \code{sommer_nonlinear_curve} specification.
The fitted object has class \code{mmes.nonlinear}, inheriting from \code{mmes}.
\code{coef()} reports ordinary fixed and natural-scale curve parameters;
\code{summary()} reports estimates and diagnostics only.
}
\examples{
set.seed(12)
dat <- data.frame(x=rep(seq(0, 5, length.out=10), 12),
                  group=factor(rep(seq_len(12), each=10)))
dat$y <- 4*exp(-0.6*dat$x) + rnorm(12, sd=0.3)[dat$group] +
  rnorm(nrow(dat), sd=0.15)
model <- nonlinear_curve("x", "exponential")
fit <- mmes(y~0+curve_term("growth"), random=~group, data=dat,
  curves=list(growth=model),
  method=mmes_method(objective="ml"),
  start=list(parameters=c(asymptote=4, rate=0.5)),
  control=mmes_control(
    convergence=list(maxIters=300),
    diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
coef(fit)
predict(fit, newdata=dat[1:3, ])
predict(fit, newdata=data.frame(x=c(1, 2)), re.form=NA)

model <- nonlinear_curve("x", random=list(
  rate=~vsm(ism(group), sigma2=0.04)))
independent <- nonlinear_curve("x", random=list(
  asymptote=~vsm(ism(group)), rate=~vsm(ism(group))))
correlated <- nonlinear_curve("x", joint=list(
  parameters=c("asymptote","rate"),
  random=~vsm(usm(.parameter),ism(group))))
# Ordinary additive effects may be added with mmes(random=~group).
# Structured rate deviations:
# nonlinear_curve("x", random=list(
#   rate=~vsm(rrm(environment, k=1), ism(group))))
}

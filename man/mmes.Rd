\name{mmes}
\alias{mmes}
\title{Mixed model equations solver}
\usage{
mmes(
  fixed,
  random = NULL,
  rcov = NULL,
  data,
  W,
  weights = NULL,
  naMethodX = "exclude",
  naMethodY = "exclude",
  naMethodRandom = "exclude",
  naMethodR = "exclude",
  contrasts = NULL,
  family = stats::gaussian(),
  vcc = NULL,
  covPar = NULL,
  parameters = list(),
  curves = list(),
  start = NULL,
  method = mmes_method(),
  control = mmes_control(),
  familyParameters = list()
)
}
\arguments{
\item{fixed}{Two-sided formula with the \strong{response(s)} and
\strong{fixed effects}, e.g. \code{y ~ Env}. Offsets are added on the
link scale. For several traits see \code{\link{stackTraits}}. The design
must be full rank (see Details).}

\item{random}{One-sided formula with the \strong{random effects}, e.g.
\code{~ Name + vsm(dsm(Env), ism(Name), Gu=Ainv)}. Terms without
\code{\link{vsm}} are identity-covariance effects. See Details for the
available covariance constructors.}

\item{rcov}{One-sided formula with the \strong{residual covariance}. The
default \code{~ units} fits one homogeneous residual variance. Structured
residuals use a single \code{vsm()} term whose last factor is normally
\code{ism(units)}, e.g. \code{~ vsm(dsm(Env), ism(units))}; complexity is
expressed by adding Kronecker factors, not by summing terms. Data need not
be sorted. In multi-trait models, replace \code{units} by a pairing key
column (e.g. \code{record} from \code{\link{stackTraits}}); without a key,
records are paired by order within trait and a warning is issued.}

\item{data}{Data frame with the model variables. Auxiliary objects
(relationship matrices, coordinates) may be found in the formula or
calling environment.}

\item{W}{Optional symmetric positive-definite observation precision (weight)
matrix, e.g. when plots are correlated. It combines with any residual
structure, including Kronecker multi-trait residuals, and with the IRLS
weights of non-Gaussian PQL fits.}

\item{weights}{Optional one-sided formula identifying independent blocks of
\code{W}, e.g. \code{~trial}. Rows from different groups must have zero
cross entries in \code{W}. It is structural metadata only and does not
generate numeric weights.}

\item{naMethodX, naMethodY, naMethodRandom, naMethodR}{Missing-value policy
for fixed predictors, responses, random-effect variables and residual
variables. The default \code{"exclude"} removes the observation from the
common mask. \code{naMethodX="include"} imputes missing covariates;
\code{naMethodY="include"} imputes missing responses and
\code{"include2"} retains multivariate rows with at least one observed
response. Imputation changes the fitted design; use with care. Supported
policies depend on the fitting route.}

\item{contrasts}{Optional named list of contrasts for factors in
\code{fixed}, as in \code{\link[stats]{lm}}.}

\item{family}{A \code{\link[stats]{family}} object or a
\code{\link{familym}} declaration. The default Gaussian identity family
fits the linear mixed model. See Details for non-Gaussian families and
estimation routes.}

\item{vcc}{Optional data frame of equality or scaling constraints among
covariance parameters (columns \code{parameter}, \code{group} and
optionally \code{scale}). Parameter names come from the \code{vcParams}
table returned with \code{mmes_control(output=list(setup=TRUE))}. See
\code{\link{vcc}}.}

\item{covPar}{Covariance parameters used as starting values, or as known
values when \code{mmes_method(task="solve")}. Usually a previous
\code{mmes} fit.}

\item{parameters}{Named list of \code{\link{parameter}()} declarations for
unknown nonlinear parameters used in \code{fixed}.}

\item{curves}{Named list of \code{\link{nonlinear_curve}()} components,
referenced in \code{fixed} with \code{\link{curve_term}("name")}.}

\item{start}{Starting values on the natural scale for nonlinear parameters
(a named vector or a list with \code{parameters}); coefficient models use
explicitly named \code{coefficients}.}

\item{method}{A \code{\link{mmes_method}()} object declaring the objective
(\code{"reml"} or \code{"ml"}), integration (\code{"auto"},
\code{"exact"}, \code{"laplace"}), approximation (\code{"pql"} or
\code{"likelihood"}) and task (\code{"fit"} or \code{"solve"}).}

\item{control}{A \code{\link{mmes_control}()} object with all numerical
settings: solver, equation formulation, convergence, output, execution
(memory, threads, GPU, precision) and diagnostics. See the
\code{sommer.mmes.options} vignette.}

\item{familyParameters}{Declarations of family nuisance parameters
(e.g. negative binomial dispersion, zero-inflation, ordinal thresholds)
as fixed or estimated.}
}
\value{
A fitted object of class \code{mmes}; its class and available
  components depend on the model and fitting route. For linear and PQL fits
  the main components are:
  \describe{
    \item{\code{b}, \code{u}, \code{bu}}{fixed effects, random effects and
      both together.}
    \item{\code{uList}, \code{uPevList}}{BLUPs and their PEVs per random
      term, with levels in rows and structure coordinates in columns.}
    \item{\code{theta}, \code{covPar}, \code{covParNative}}{estimated
      covariance matrices, parameters in reported coordinates, and a table
      of parameters on the native model scale with standard errors and
      z-ratios.}
    \item{\code{theta_se}, \code{InfMat}}{covariance-parameter uncertainty
      and information matrix; \code{NULL} when
      \code{covarianceInformation="none"}.}
    \item{\code{Ci}, \code{Ci_11}, \code{Cchol}}{inverse and Cholesky factor
      of the coefficient matrix.}
    \item{\code{llik}, \code{AIC}, \code{BIC}, \code{convergence},
      \code{monitor}}{likelihood, information criteria, convergence flag and
      the iteration history.}
    \item{\code{data}, \code{obsInfo}, \code{y}, \code{Dtable}}{retained
      data, the missing-data audit, the response and the factor table used
      by \code{predict()}.}
    \item{\code{partitions}, \code{partitionsX}}{start and end positions of
      each random and fixed term.}
    \item{\code{toBoundary}, \code{remlOmissions}}{components forced to the
      boundary and structurally unestimable components omitted from REML.}
    \item{\code{engineDiagnostics}}{path-dependent diagnostics
      (compression, PCG, profile timings).}
    \item{\code{control}}{validated controls used in the fit.}
  }
  Use the extractors and methods \code{summary()}, \code{coef()},
  \code{fitted()}, \code{residuals()}, \code{predict()}, \code{anova()},
  \code{vcov()}, \code{plot()}, \code{\link{randef}},
  \code{\link{covparams_mmes}}, \code{\link{covparams_mmes_se}},
  \code{\link{postPEV}}, \code{\link{postVarU}}, \code{\link{wald_mmes}},
  \code{\link{boot_mmes}}, \code{\link{r2}} and \code{\link{vcDelta}}.
}
\description{
\code{mmes()} is the main model-fitting function of sommer. It fits linear,
generalized and nonlinear mixed models with arbitrary Kronecker-structured
random and residual covariance, known relationship matrices, multiple
responses, spatial structures and low-rank (factor analytic) models.
Variance components are estimated by REML or ML using sparse or dense
average-information, EM, PX-EM and related solvers implemented in C++.
The same call can also solve the mixed model equations for known variance
components (genetic evaluation).
}
\details{
\strong{Specification layers.} A model is built from four independent
pieces:
\enumerate{
  \item \emph{What is modeled}: the formulas \code{fixed}, \code{random}
    and \code{rcov}, plus \code{family}.
  \item \emph{Which covariance structures}: \code{\link{vsm}()} terms inside
    \code{random} and \code{rcov} (see below).
  \item \emph{What is estimated}: the statistical objective, integration and
    task, declared with \code{\link{mmes_method}()}.
  \item \emph{How it is computed}: solver, formulation, convergence,
    output, threading and diagnostics, declared with
    \code{\link{mmes_control}()}.
}

\strong{Covariance structures with vsm().} Every structured random or
residual term is written as one \code{\link{vsm}()} call:
\deqn{G = \sigma^2 (K_1 \otimes K_2 \otimes \cdots \otimes K_m) \otimes A,}
where \eqn{\sigma^2} is the single variance scale owned by \code{vsm()},
\eqn{K_1,\ldots,K_m} are dimensionless covariance shapes built by the
constructors listed below, and \eqn{A} is an optional known relationship
matrix among the levels of the last term. The \strong{last} argument of
\code{vsm()} supplies the incidence matrix of the main effect; all previous
arguments are covariance factors. There is no limit on the number of
factors. A random term not wrapped in \code{vsm()} is an
identity-covariance effect, so \code{random = ~ Name} equals
\code{random = ~ vsm(ism(Name))}.

Known relationships are given through \code{Gu}, the \strong{inverse} of
\eqn{A}, with row and column names equal to the levels of the last term and
\code{attr(Gu, "inverse") = TRUE}; see \code{\link{vsm}} for the accepted
formats. Matrices can be computed with \code{\link{A.mat}} and
\code{\link{H.mat}}.

Covariance constructors available:
\tabular{ll}{
  \code{\link{ism}} \tab identity \cr
  \code{\link{dsm}}, \code{\link{atm}} \tab diagonal (heterogeneous
    variances), diagonal over selected levels \cr
  \code{\link{usm}} \tab unstructured \cr
  \code{\link{csm}}, \code{\link{corgm}} \tab compound symmetry and
    correlation structures \cr
  \code{\link{ar1m}}, \code{\link{ar2m}}, \code{\link{ar3m}} \tab
    autoregressive \cr
  \code{\link{mam}} (\code{ma1m}, \code{ma2m}) \tab moving average \cr
  \code{\link{toeplitzm}} \tab Toeplitz correlation \cr
  \code{\link{antem}} \tab antedependence \cr
  \code{\link{fam}} \tab factor analytic (see the factor analytic vignette)
    \cr
  \code{\link{rrm}} \tab reduced rank, specific variances fixed at zero.
    Use \code{fam(..., variance="homogeneous")} for a common estimated
    specific variance \cr
  \code{\link{maternm}} \tab Matern spatial covariance \cr
  \code{\link{metricm}} \tab distance-based exponential, power, gaussian,
    spherical and circular correlations (isotropic or anisotropic) \cr
  \code{\link{sar}}, \code{\link{car}} \tab simultaneous and conditional
    autoregressive spatial models \cr
  \code{\link{ownm}} \tab user-defined or fixed covariance \cr
  \code{\link{dsumm}} \tab direct sum: separate structure per group \cr
}
Design-matrix utilities such as \code{\link{spl2Dc}},
\code{enhancer::overlay}, \code{enhancer::leg} and \code{enhancer::redmm}
can be combined with \code{vsm()}. Variance components can be constrained
with \code{\link{vcc}}, \code{\link{fixm}} or \code{\link{covm}}.

Examples of \code{random} and \code{rcov} specifications:
\tabular{ll}{
  \code{~ Name + Env:Name} \tab compound symmetry (main plus interaction) \cr
  \code{~ vsm(dsm(Env), ism(Name))} \tab heterogeneous genetic variance by
    environment \cr
  \code{~ vsm(usm(Env), ism(Name), Gu=Ainv)} \tab unstructured GxE with
    known relationships \cr
  \code{~ vsm(fam(Env, k=2), ism(Name), Gu=Ainv)} \tab factor analytic GxE
    \cr
  \code{~ vsm(ar1m(Row), ar1m(Col), ism(Plot))} \tab separable spatial AR1 \cr
  \code{rcov = ~ vsm(dsm(Env), ism(units))} \tab heterogeneous residual
    variances \cr
  \code{rcov = ~ vsm(usm(trait), ism(record))} \tab multi-trait residual
    with explicit pairing key \cr
  \code{rcov = ~ dsumm(vsm(ar1m(Row), ar1m(Col), ism(units)), by=Trial)}
    \tab section-specific spatial residual \cr
}

\strong{Non-Gaussian responses.} Non-Gaussian mixed models use penalized
quasi-likelihood (PQL) by default: IRLS working responses are fitted with
the Gaussian REML machinery. Because the working likelihood is not a GLMM
likelihood, \code{AIC}, \code{BIC} and \code{llik} are \code{NA} and
likelihood-ratio \code{anova()} is refused for PQL fits; use deviance and
Wald tests (\code{\link{wald_mmes}}). Use
\code{method=mmes_method(objective="ml", integration="laplace")} for the
native Laplace likelihood. Supported families include binomial (logit,
probit, cloglog), Poisson, Gamma, inverse Gaussian, Gaussian with other
links, quasi families, \code{\link{nbinom2}}, \code{\link{zipoisson}},
\code{\link{ordinal}}, \code{\link{multinomial}}, grouped binomials via
\code{\link{binm}}, and mixtures through \code{\link{familym}}. Fixed-only
generalized models require \code{objective="ml"}.

\strong{Nonlinear models.} Unknown nonlinear parameters declared in
\code{parameters} with \code{\link{parameter}()}, or native curves in
\code{curves} declared with \code{\link{nonlinear_curve}()} and used in the
formula through \code{\link{curve_term}()}, require
\code{mmes_method(objective="ml")}. Fits are post-processed with
\code{\link{profile_nonlinear}} and \code{\link{post_nonlinear}}.

\strong{Known variance components.} With
\code{method=mmes_method(task="solve")} and \code{covPar} set to a previous
fit, the mixed model equations are solved without estimating variances
(genetic evaluation, see the genetic evaluation vignette).

\strong{Missing data.} All model components are evaluated before a common
observation mask is applied, so response, fixed, random and residual
designs always use the same rows. The \code{naMethod*} arguments set the
policy per component and the audit is returned in \code{obsInfo}.

\strong{Fixed-effect rank.} The fixed design must have full column rank.
Aliased or all-zero columns are reported before fitting. Request
\code{mmes_control(representation=list(fixedRank="drop"))} to drop aliased
coefficients (returned as zeros and flagged in \code{summary(fit)$betas}),
or use \code{droplevels()} and \code{interaction(..., drop=TRUE)} to remove
unobserved cells.

\strong{Rotation.} \code{vsm(..., Gu=Gu, rotation=TRUE)} requests
Lee--van der Werf rotation for one relationship term. It requires a
balanced Gaussian layout and a rotation-invariant residual structure.

\strong{Where to learn more.} Vignettes (\code{browseVignettes("sommer")}):
\itemize{
  \item \code{vignette("sommer.models")}: model specification.
  \item \code{vignette("sommer.covariance.structures")}: all covariance
    structures and \code{vsm()} usage.
  \item \code{vignette("sommer.gxe")}: genotype-by-environment models.
  \item \code{vignette("sommer.fa.and.rr.models")}: factor analytic and
    reduced rank models.
  \item \code{vignette("sommer.spatial")}: spatial models.
  \item \code{vignette("sommer.qg")}: quantitative genetics.
  \item \code{vignette("sommer.genetic.evaluation")}: known variance
    components.
  \item \code{vignette("sommer.glmm")} and
    \code{vignette("sommer.glmm.theory")}: generalized models.
  \item \code{vignette("sommer.mmes.options")}: \code{mmes_control()} and
    \code{mmes_method()} options.
  \item \code{vignette("sommer.vs.lme4")}: translating lme4 models.
  \item \code{vignette("sommer.theory")} and
    \code{vignette("sommer.machinery")}: parameterization and engine.
}
}
\examples{
dat <- data.frame(x=seq_len(12), y=1+seq_len(12)/2+sin(seq_len(12))/10)
fit <- mmes(y~x, data=dat, control=mmes_control(
  diagnostics=list(verbose=FALSE, dateWarning=FALSE)))
coef(fit)

\donttest{
data(DT_example, package="enhancer")
DT <- DT_example

## compound symmetry: genotype main effect + genotype-by-environment
ans1 <- mmes(Yield ~ Env, random = ~ Name + Env:Name,
             rcov = ~ units, data = DT)
summary(ans1)

## known variance components: solve the MME only
sol1 <- mmes(Yield ~ Env, random = ~ Name + Env:Name, rcov = ~ units,
             data = DT, covPar = ans1,
             method = mmes_method(task = "solve"))
head(sol1$uList[[1]])

## heterogeneous genetic and residual variances by environment
DT <- DT[order(DT$Env), ]
ans2 <- mmes(Yield ~ Env,
             random = ~ Name + vsm(dsm(Env), ism(Name)),
             rcov = ~ dsumm(vsm(ism(units)), by = Env), data = DT)
summary(ans2)

## unstructured genotype-by-environment covariance
ans3 <- mmes(Yield ~ Env, random = ~ vsm(usm(Env), ism(Name)),
             rcov = ~ dsumm(vsm(ism(units)), by = Env), data = DT)
summary(ans3)

## binomial GLMM by PQL (default) and Laplace likelihood
set.seed(1)
g <- data.frame(id = factor(rep(1:30, each = 6)), x = rnorm(180))
g$y <- rbinom(180, 1, plogis(0.5*g$x + rep(rnorm(30), each = 6)))
pql <- mmes(y ~ x, random = ~ id, data = g, family = binomial())
lap <- mmes(y ~ x, random = ~ id, data = g, family = binomial(),
            method = mmes_method(objective = "ml", integration = "laplace"))
summary(lap)
}
}
\references{
Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits
using the R package sommer. PLoS ONE 2016, 11(6):
doi:10.1371/journal.pone.0156744

Jensen J, Mantysaari EA, Madsen P, Thompson R. 1997. Residual maximum
likelihood estimation of (co)variance components in multivariate mixed
linear models using average information. Journal of the Indian Society of
Agricultural Statistics 49:215-236.

Gilmour AR, Thompson R, Cullis BR. 1995. Average information REML: an
efficient algorithm for variance parameter estimation in linear mixed
models. Biometrics 51(4):1440-1450.

Sanderson C, Curtin R. 2025. Armadillo: an efficient framework for
numerical linear algebra. arXiv:2502.03000.
}
\seealso{
\code{\link{vsm}} and the constructors \code{\link{ism}},
  \code{\link{dsm}}, \code{\link{usm}}, \code{\link{csm}},
  \code{\link{ar1m}}, \code{\link{fam}}, \code{\link{rrm}},
  \code{\link{maternm}}, \code{\link{ownm}}; \code{\link{mmes_method}},
  \code{\link{mmes_control}}; \code{\link{summary.mmes}},
  \code{\link{predict.mmes}}, \code{\link{anova.mmes}},
  \code{\link{wald_mmes}}; \code{\link{familym}}, \code{\link{parameter}},
  \code{\link{nonlinear_curve}}; \code{\link{A.mat}}, \code{\link{H.mat}},
  \code{\link{stackTraits}}, \code{\link{vcc}}
}
\author{
Giovanny Covarrubias-Pazaran
}

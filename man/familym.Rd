\name{familym}
\alias{familym}
\title{Different distribution families for the traits of a multi-trait model}
\description{
Assigns a \code{stats::family()} object to each trait of a long-format
multi-trait \code{\link{mmes}} model, e.g. a binary disease score together
with a Gaussian yield. PQL supports its existing working-model covariance
structures; native Laplace supports a joint latent-effect likelihood with
conditionally independent observations.
}
\usage{
familym(by, ..., trials = NULL)
}
\arguments{
  \item{by}{Name (character string) of the column of \code{data} that
  identifies the traits, typically \code{"trait"} as created by
  \code{\link{stackTraits}}.}
  \item{\dots}{Family objects named by the levels of the \code{by} column,
  e.g. \code{Disease = binomial(), Yield = gaussian()}. Every level present in
  the data needs a family.}
  \item{trials}{Optional trials-column name for native Laplace fits. Binomial
  responses are then integer success counts, with positive integer trial
  totals in this column. All non-binomial rows must have a trial value of one.
  Omit this argument for binary responses or responses made with \code{binm()}.
  Do not declare trials both here and in a response attribute/matrix.}
}
\details{
The object is passed as the \code{family} argument of \code{mmes}. Each PQL
iteration builds the working response and working weights of every row with
the family of its trait, and fits a weighted Gaussian multi-trait model.

The residual variance of binomial, Poisson and negative binomial traits is not
a free parameter (their dispersion is one), while that of Gaussian and other
families is estimated. The residual structure must therefore allow individual
trait variances to be fixed: use
\code{rcov = ~ vsm(corgm(trait, variance = "heterogeneous"), ism(record))}
(trait correlations and heterogeneous variances, see \code{\link{corgm}}) or
\code{rcov = ~ vsm(dsm(trait), ism(units))} (uncorrelated residuals).
The first level of the trait factor must be a fixed-dispersion trait,
because its variance is the \code{vsm} scale; reorder the levels
otherwise. \code{mmes} then fixes the residual variance of every
fixed-dispersion trait to one and estimates the rest.

As for single-family PQL fits, the likelihood, \code{AIC} and \code{BIC} are
\code{NA}. The fitted object has a named \code{dispersion} vector (one per
trait: one for fixed-dispersion traits, the REML residual variance for the
others), and \code{fitted()} and \code{residuals()} apply the inverse link and
variance function of each trait to its rows. Random-effect variances are on
the link scale of each trait.

For native Laplace, select \code{task=list(objective="ml",
glmm=list(method="laplace"))}. Supported families are binomial with logit,
probit or complementary log-log link; Poisson-log;
\code{\link{zipoisson}} with an intercept-only or fixed-effects formula
logit structural-zero component; \code{nbinom2()}-log; Gaussian-identity;
and Gamma-log. Single-family ordered categorical models are available through
\code{\link{ordinal}} with cumulative logit or probit links. Ordinal
\code{familym()} models remain explicitly unsupported.
Zero-inflated Poisson is not a hurdle model and its zero component has no
random effects. Observations are conditionally independent
given the shared random effects. A correlated \code{usm(trait)} random
covariance can connect traits, but structured residual covariance is not
supported on this route: omit \code{rcov} or use \code{~units}. There is no
requirement to place a fixed-dispersion trait first.

NB2 shape and Gaussian/Gamma dispersion are estimated separately for each
declared trait by default; zero-inflated Poisson structural-zero probability
is also estimated by default. Configure starting/fixed values under
\code{task$glmm$familyParameters}, for example
\code{list(count=list(theta=2, estimate=FALSE),
zeroCount=list(zeroInflation=0.2, estimate=TRUE),
continuous=list(dispersion=0.7, estimate=TRUE))}.
Every declared family must retain observations. Missing trait/trials values
participate in the fixed-covariate missing-data policy; unknown non-missing
trait labels error. The \code{trials} column argument is currently
Laplace-only; existing PQL trial responses use \code{binm()}.

Native fits report an approximate marginal log likelihood, AIC/BIC with
estimated shape, dispersion or zero-inflation parameters counted, a
\code{familyParameters} table for scalar family parameters,
\code{zeroCoefficients} for formula-based zero components, per-trait
\code{dispersion} and observation-family indices. Zero-only missing covariates
affect only their declared ZIP trait's rows. Fitted
values and residuals apply each row's fitted family. New-data response-scale
predictions require the declared trait column. Use \code{re.form=NA} for
zero-random-effect means, \code{re.form=NULL} for conditional predictions at
known random-effect levels, or \code{re.form="population"} for supported
population means; see \code{\link{predict.mmes}} for restrictions.
}
\value{An object of class \code{"sommer_familym"}.}
\seealso{\code{\link{mmes}}, \code{\link{corgm}}, \code{\link{stackTraits}},
\code{\link{covmatrix_mmes}}}
\examples{
set.seed(12)
ng <- 100; nr <- 4
G <- matrix(c(1, 0.4, 0.4, 0.8), 2)
u <- matrix(rnorm(ng * 2), ng) \%*\% chol(G)
d <- data.frame(id = factor(rep(paste0("g", 1:ng), each = nr)))
d$disease <- rbinom(nrow(d), 1, plogis(-0.5 + u[as.integer(d$id), 1]))
d$yield <- 10 + u[as.integer(d$id), 2] + rnorm(nrow(d))
L <- stackTraits(d, traits = c("disease", "yield"))

\donttest{
fit <- mmes(value ~ trait,
  random = ~ vsm(usm(trait), ism(id)),
  rcov = ~ vsm(corgm(trait, variance = "heterogeneous"), ism(record)),
  family = familym("trait", disease = binomial(), yield = gaussian()),
  data = L,
  control = mmes_control(diagnostics = list(verbose = FALSE)))
fit$dispersion
covmatrix_mmes(fit, 1)$correlation   # genetic correlation (link scale)
head(fitted(fit))                    # probabilities for disease rows
}


}
\keyword{models}

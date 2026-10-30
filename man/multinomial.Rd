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
uncertainty is deferred until after fitting. PQL, \code{familym()}, supplied
\code{W}, weights formulas, \code{vcc}, structured residual covariance,
and non-default Gaussian solver/acceleration controls are unsupported.
Random-effect relationships and covariance structures use the existing
\code{mmes_control()} and covariance-descriptor machinery. No new package
dependency is required.

\code{fitted()} and \code{predict()} return a matrix of category probabilities
with the original response levels as columns; \code{type="link"} returns the
\eqn{K-1} logits. Conditional prediction uses fitted modes for known random
levels, and \code{re.form=NA} sets random effects to zero. Population-marginal
prediction is not yet supported. \code{simulate()} returns factors preserving
the original category labels and can draw new population random effects with
\code{re.form="population"}.

Post-fit \code{summary()} and \code{vcov()} calculate fixed-effect uncertainty
from the numerical Hessian of the marginal Laplace objective. \code{postPEV()}
calculates random-effect posterior variances after fitting.
\code{predict(se.fit=TRUE)} uses the joint fixed/random conditional curvature;
its probability-scale standard errors use a delta-method transformation.
These prediction standard errors condition on the fitted covariance
parameters. Boundary, nonconverged, or non-positive-definite fits do not
receive success-shaped uncertainty. \code{boot_mmes()} supports nominal
parameter-recovery and nested fixed-effect refits by simulating population
random effects and refitting the joint softmax likelihood. Bootstrap draws
are reported on the fixed-coefficient and native covariance-coordinate
scales. Failed or nonconverged replicates remain in diagnostics, and
percentile intervals use successful refits only. This descriptive bootstrap
has not yet been calibrated as a replacement for Wald intervals.

Likelihood, AIC and BIC are available for the fitted nominal model, and
nested fixed-effect likelihood comparisons are supported when their response,
random structure and covariance specification match. Statistical calibration
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
summary(fit)$betas                # post-fit fixed-effect uncertainty
simulate(fit, nsim=2, seed=10)
}
}
\keyword{models}

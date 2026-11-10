\name{parameter}
\alias{parameter}
\title{Declare a nonlinear predictor parameter regression}
\usage{
parameter(fixed = ~1, random = NULL, link = "identity")
}
\arguments{
\item{fixed}{A one-sided fixed coefficient regression, for example \code{~Env}.}
\item{random}{A one-sided covariance formula using the common \code{vsm} machinery.}
\item{link}{\code{"identity"} or \code{"log"}, mapping the natural parameter to its regression predictor.}
}
\description{
Declares a named parameter in a model formula expression. Parameter random
effects are added before applying the inverse parameter link. A log-linked
rate is therefore positive, independently of the response family's link.
}
\value{A \code{sommer_parameter} declaration, not a fitted model.}
\details{
Nonlinear specification models require maximum likelihood. Gaussian identity
models with only additive random effects use exact integration; parameter
random effects and non-Gaussian random effects use conditional-mode Laplace
integration with numerical scores and Hessians. This bounded dense reference
implementation supports Gaussian, Bernoulli, Poisson, and Gamma-log likelihoods,
independent residuals (\code{rcov=NULL} or \code{rcov=~units}),
up to 60 mean coefficients, and up to 512 random
coefficients for Laplace integration. It never adds artificial random effects
or regularizes a covariance.
Outer optimization uses \code{nlminb}; conditional modes use numerical BFGS.
Numerical conditional scores must meet the requested mode tolerance;
step tolerances detect stalled mode refinement. Non-default solver/update
modules, weights, custom contrasts, covariance constraints, nonempty
\code{familyParameters}, and post-inference controls are explicitly
unsupported. Automatic or central outer gradients are numerical; an explicit
analytic-gradient request is rejected. Iteration limits, convergence
tolerances, and population quadrature controls are supported.

Expressions support arithmetic, powers, \code{exp}, \code{log}, \code{tanh},
\code{sqrt}, \code{sin}, and \code{cos}; data inputs must be finite numeric
columns. Ordinary additive formula terms retain standard R formula semantics.
Unused declarations are errors. The \code{naMethodX}, \code{naMethodY},
\code{naMethodRandom}, and \code{naMethodR} policies accept \code{"exclude"}
or \code{"fail"}. A shared complete-observation mask is applied across the
response, expression inputs, regressions, and random grouping columns.
Fitted values and residuals use the retained rows, recorded in \code{obsInfo}.

Natural starts such as \code{start=c(a=3,k=1)} initialize parameter intercepts;
other parameter regression coefficients start at zero. Alternatively use
\code{start=list(parameters=c(a=3,k=1))}. Regressions without an intercept
require \code{start=list(coefficients=...)} containing every mean coefficient.
Those names are \code{fixed:(Intercept)}, \code{a:(Intercept)}, or
\code{k:EnvB}, and values are on coefficient/link scales.

Prediction distinguishes conditional, zero-random-effect, and population
targets, on response or link scales. Population prediction uses adaptive
Gaussian quadrature and supports at most two independent random predictor
coordinates. Unknown random levels are rejected for conditional and population
prediction. Coefficient covariance is computed only when \code{vcov()} is
requested and requires positive-definite observed numerical curvature.
Fitting and summary do not compute outer parameter uncertainty.
Simulation is explicitly unsupported.
}
\examples{
parameters <- list(a=parameter(), k=parameter(fixed=~1, link="log"))
parameters$k
}
\seealso{\code{\link{graph_term}}, \code{\link{vsm}}}

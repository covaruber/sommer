\name{vcc}
\alias{vcc}
\alias{vcParams}
\title{Equality and scaling constraints between variance parameters}
\description{
The \code{vcc} argument of \code{\link{mmes}} constrains groups of covariance
parameters to be equal or to keep fixed ratios.
Parameters are identified through the \code{vcParams} table.
}
\section{The vcParams table}{
\code{mmes(..., control = mmes_control(output = list(setup = TRUE)))$vcParams}
(before fitting) and
\code{fit$vcParams} (after fitting) list every covariance parameter with
columns \code{index}, \code{term}, \code{parameter}, \code{kind},
\code{start} (natural scale), \code{free}, \code{structure} and \code{position}
(and \code{group} when constraints were used). \code{kind} is \code{"scale"}
for the variance scale (\code{sigma2}) owned by each \code{vsm()} term, or the
transform used to estimate the parameter: \code{"exp"} (positive quantities
such as ranges and variance ratios), \code{"identity"} (e.g. Cholesky or
factor-analytic loadings), \code{"tanh"} or \code{"bounded_logit"}
(correlation-type parameters).
}
\section{Specifying constraints}{
\code{vcc} is a data frame with columns
\describe{
  \item{parameter}{the \code{index} of the parameter, or its name
  \code{"term:parameter"} (e.g. \code{"vsm(ism(B)):sigma2"}); the
  \code{parameter} name alone is accepted when unique.}
  \item{group}{parameters with the same group value are constrained together.}
  \item{scale}{optional, default 1: members divided by their respective
  \code{scale} values are equal on the natural scale. Relative to the first
  member, the ratio is \code{scale / scale[1]}. Scaling is allowed for
  \code{"scale"}, \code{"exp"} and \code{"identity"} parameters; correlation-type
  parameters can only be equated.}
}
All members of a group must have the same \code{kind}: variance scales can be
constrained with other variance scales, and other parameters with parameters
estimated through the same transform. If any member of a group is fixed, the
whole group is fixed at the value implied by that member. Starting values of
the members are made consistent with the first (or fixed) member.
Multiple fixed members must imply the same value. Correlations with bounded
logit transforms must have identical bounds to be equated.
}
\details{
Parameters are estimated on transformed scales that guarantee valid
covariance matrices (log for variances, atanh or bounded logits for
correlations). A constraint is enforced exactly when it is affine on that
scale, which covers equalities of any parameter and ratios of positive
parameters (a log-scale offset). Relations such as
\eqn{\sigma^2_3=\sigma^2_1+\sigma^2_2} are not affine on the log scale and are not supported.

Internally the free working parameters are written as \eqn{w=o+T\phi}. Each
average-information update is projected onto this set in the information
metric, which equals the constrained Newton step, and then repaired exactly;
both the Henderson and the direct engines support constraints. The covariance
of the estimates is \eqn{T(T^\top\mathcal{I}T)^{-1}T^\top}, so constrained
parameters share their standard errors.
}
\section{Native Laplace GLMMs}{
Native Laplace ML supports these constraints for random covariance parameters,
including ordinal, nominal and mixed \code{familym()} responses. Use the same
family and \code{method$integration="laplace"} when requesting
\code{output$setup=TRUE}: the returned \code{vcParams} table contains random
covariance parameters, plus actual Gaussian observation-residual parameters
when structured Gaussian blocks are requested. The latter are reported for
inspection and uncertainty but are not accepted as \code{vcc} members.
The synthetic Gaussian preparation residual scale is always excluded.
Parameter indices from a Gaussian setup are not interchangeable with this table.
Family dispersion, NB2 shape, thresholds and zero-inflation coefficients are not
part of \code{vcc}.

The likelihood is optimized in the independent coordinates \eqn{\phi}; native
scores are pulled back with \eqn{T^\top}. Parameter counts for \code{logLik()},
AIC and BIC count only independent estimated coordinates. Numerical observed
information, covariance standard errors, prediction uncertainty and PEV remain
post-fit calculations. Their Jacobians propagate uncertainty through the same
constraint map and the natural-scale transformations; no Hessian is computed
during fitting.

Likelihood comparisons through \code{anova()} still require nested fixed
effects with the same covariance constraints. Comparing different variance
constraint spaces is not supported by this fixed-effects comparison.
}
\seealso{\code{\link{mmes}}, \code{\link{vsm}}, \code{\link{covparams_mmes}}}
\examples{
data(DT_yatesoats, package = "enhancer")
DT <- DT_yatesoats
p <- mmes(Y ~ V * N,
  random = ~ B + B:MP,
  rcov = ~ units,
  data = DT,
  control = mmes_control(output = list(setup = TRUE)))
p$vcParams
# block and whole-plot variances equal
m1 <- mmes(Y ~ V * N,
  random = ~ B + B:MP,
  rcov = ~ units,
  data = DT,
  vcc = data.frame(parameter = c(1, 2), group = 1),
  control = mmes_control(diagnostics = list(verbose = FALSE)))
# block variance twice the whole-plot variance
m2 <- mmes(Y ~ V * N,
  random = ~ B + B:MP,
  rcov = ~ units,
  data = DT,
  vcc = data.frame(parameter = c("vsm(ism(B:MP)):sigma2", "vsm(ism(B)):sigma2"),
                            group = 1, scale = c(1, 2)),
  control = mmes_control(diagnostics = list(verbose = FALSE)))
unlist(m2$covPar)
}
\keyword{models}

\name{profile_nonlinear}
\alias{profile_nonlinear}
\title{Profile traditional nonlinear mixed-model operations}
\description{
Reconstructs operation-local native contexts and measures separate nonlinear
evaluation, derivative, conditional-curvature, marginal-likelihood and prediction
operations. Does not refit or mutate the supplied model.
}
\usage{
profile_nonlinear(object, replications=3L, control=NULL, inference=FALSE)
}
\arguments{
\item{object}{A fitted traditional nonlinear \code{mmes} model.}
\item{replications}{Positive integer number of timed repetitions after one warmup.}
\item{control}{Optional \code{mmes_control()} providing the profiling memory
budget and, when requested, post-fit inference settings. Other fitted
numerical choices remain unchanged.}
\item{inference}{Whether to explicitly profile post-fit uncertainty as well.
Defaults to \code{FALSE}; profiling never implicitly computes PEV.}
}
\details{
Use \code{mmes_control(diagnostics=list(profile="force"))} at fitting time
to retain preparation, outer-optimization and result-construction elapsed
seconds. Marginal evaluation time is included in outer optimization and
must not be added again. Default fitting retains no timing profile.

The profiling table measures curve means, curve Jacobians, exact conditional
curvature at the saved latent mode, the marginal objective (including mode
search), the marginal objective with its analytic score, and conditional
prediction. These operation timings are inclusive, overlapping and not
additive. The difference between objective and score timings can indicate
derivative cost, but timer resolution and run-to-run variation limit that
comparison. Context reconstruction is timed separately. Whole-optimizer
timing and iteration/evaluation counts come from the fit, not a second fit.
No automatic inference, stochastic starts or backpropagation is introduced.

Reported memory consists of retained R object sizes and the actual
workspace admission policy, not measured native peak resident memory.
Before nonlinear Hessian assembly a sparse-support graph bounds dense fill
separately within each connected random component, including cross-category
curvature edges. This bound scales linearly for independent groups but can
be conservative for a large connected sparse relationship graph.
The estimated workspace is
\deqn{8(180n+80q+12\,nnz(Z)+8\,nnz(Q)+8\,nnz(C))
       +64\sum_j q_j^2}
bytes, with predictor rows \eqn{n}, latent coordinates \eqn{q}, prior
precision \eqn{Q}, observation curvature \eqn{C}, and component sizes
\eqn{q_j}. Admission occurs before sparse Hessian products/factorization.
For random curves the support includes both additive and curve incidence,
even when a curve derivative is zero. This is a conservative operation-local
bound, not an operating-system cap or an exact symbolic-fill estimate.
Existing score-panel, reduced-rank, prediction, quadrature and post-fit dense
workspace checks also apply. No full observation or random design is
densified to evaluate the curve.
}
\value{
A serializable list with \code{table}, context reconstruction time,
fit-time profile (if requested when fitting), dimensions, cumulative native
work counters, optimizer/initialization diagnostics, memory information and
R/platform/BLAS metadata. Timings are seconds; \code{fitBytes} and
\code{sparseDesignBytes} are bytes. Very short operations can time as zero.
}
\seealso{\code{\link{nonlinear_curve}}, \code{\link{post_nonlinear}},
\code{\link{mmes_control}}}
\examples{
\dontrun{
fit <- mmes(y~0+graph_term("growth"), random=~group, data=dat,
  predictors=list(growth=nonlinear_curve("x")),
  method=mmes_method(objective="ml"),
  control=mmes_control(
    diagnostics=list(profile="force")))
profile_nonlinear(fit)
}
}

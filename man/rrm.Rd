\name{rrm}
\alias{rrm}
\title{Reduced-rank covariance structure}
\description{
Reduced-rank covariance structure. The returned dimensionless covariance shape is intended for use inside \code{vsm}, which supplies the single overall variance scale.
}
\usage{
rrm(x, k = 1L, loadings = NULL, fixed = NULL)
}
\arguments{
  \item{x}{Variable defining q covariance levels.}
  \item{k}{Reduced rank, satisfying 1 <= k < q.}
  \item{loadings}{Optional finite q by k starting loading matrix. The leading
k by k block is made lower triangular and must have positive diagonal entries.
Loadings are divided by their first diagonal entry to remove the internal scale.}
  \item{fixed}{Logical vector of length \eqn{qk-k(k-1)/2}, controlling the
triangular loading entries in column-major order. The first entry is always
fixed: the reference loading is one, and \code{vsm()} estimates the overall
variance scale. Other fixed entries constrain relative loadings.}
}
\details{
The reduced-rank structure uses
\deqn{K=\Lambda\Lambda^{\mathsf T},\qquad \Lambda_{11}=1,}
and \code{vsm()} supplies the overall variance scale:
\deqn{G=\sigma^2 K.}
All specific variances are fixed at exactly zero; no diagonal remainder or
numerical nugget is added. This is the same covariance family as ASReml's
default \code{rr(x, k)} model. Individual loadings can differ by sign or
rotation, so compare reconstructed covariance matrices when comparing fits.

The leading \eqn{k\times k} loading block is lower triangular for rotational
identification and its diagonal entries are positive. There are
\eqn{qk-k(k-1)/2-1} free shape parameters plus the overall scale.
Because \eqn{k<q}, the random-effect covariance is positive semidefinite
and singular by design.

The Henderson engine fits this structure through latent factors rather than
inverting the singular environment covariance. The direct-inversion formulation
(\code{control$formulation="direct"}) also supports it, provided the total observation
covariance is positive definite. \code{rrm()} is intended for random effects,
not residual covariance structures.

Native latent-factor Henderson fits support \code{control$solve$method="ldlt"} and
\code{control$solve$method="cholmod"} (including automatic selection). They use
loading-aware likelihood scores and average-information updates in the
original loading/variance coordinates. The loading derivatives account for
the parameter-dependent latent design; the singular environment covariance
is never inverted. \code{control$proposal$stepWeight} damps covariance-parameter proposals,
while \code{control$update$blendWeight} blends EM/GEM or scalar PX-GEM with
average-information update directions, selected through \code{control$update$recovery},
as in ordinary \code{mmes()} fitting. Parameter standard errors use the unblended
average-information matrix at the final parameter values.
Prediction errors are mapped back from the latent coefficient equations
to the original random effects.

The likelihood trace terms are computed exactly by sharing random-design
and residual-block contractions across covariance parameters. Bounded panels
of Henderson solves avoid constructing a dense record-by-record covariance
or projection matrix. Random trace solves use the loading span (rank \code{k}),
rather than every environment direction. Reducing the coefficient dimension does not guarantee
a faster fit for every dataset.

Section-specific spatial residuals supplied by \code{dsumm()} are evaluated
on their observed residual blocks, including grids with missing cells and
unlisted sections with independent residuals. Residual covariance derivatives
and trace contractions never expand the full section-by-grid covariance.
\code{engineDiagnostics$residualEvaluation} reports \code{"observed-blocks"}
and \code{engineDiagnostics$largestResidualBlock} records the largest observed
block dimension used by the native RR AI evaluator.

Likelihood and AI calculations share the selected residual precision operator
and the live coefficient factorization, including grid residuals with
missing-cell absorption and eligible block-Schur solves. No separate
coefficient factorization or dense observed-block precision is built for AI.
Random trace panels use the latent Henderson identity to cancel relationship
covariance actions. Grid residual traces use block-local projection panels
rather than repeated full-observation precision applications.
\code{engineDiagnostics$residualEngine} identifies the
shared residual engine, and \code{sharedConditionalOperators} confirms its
use. \code{aiElapsedSeconds} gives elapsed timings for the final AI
evaluation's setup, random traces, residual traces, and information matrix.

The positive-specific-variance Woodbury Schur formulation used for
\code{fam()} is unnecessary for exact RR: the specific-effect blocks are
absent, not inverted at zero. The smaller latent coefficient system still
uses the ordinary LDLT/CHOLMOD and eligible block-Schur dispatch.
With \code{control$diagnostics$verbose=TRUE}, the native fit reports its RR path, selected
coefficient solver, and then the progress table. Each accepted iteration
reports likelihood, wall-clock time, cumulative CPU seconds, the scheduled
\code{EMweight}, accepted-step \code{trustAlpha} and \code{lineAlpha},
and score norm. The step scales describe the update that reached the current
row (both are one at initialization); \code{lineAlpha} includes
\code{control$proposal$stepWeight} and any backtracking. Boundary-restraint and \code{notPD}
columns are omitted because the latent optimizer does not use those
ordinary covariance-update diagnostics. Numeric weights, step scales and
CPU times are retained in \code{engineDiagnostics$iterationProgress}.

For native Henderson fits with weights \code{W}, let
\eqn{U=\mathrm{chol}(W)} in observation order. The effective residual
precision is \eqn{U R_0^{-1} U^{\mathsf T}}, where \eqn{R_0} is the
unweighted residual covariance. The score and AI use derivatives of the
corresponding covariance \eqn{U^{-\mathsf T}R_0U^{-1}}.

The native path does not currently support PCG, \code{vcc} constraints,
covariance-update acceleration, \code{control$task$mode="solve"}, or
\code{control$representation$factorScores} augmentation. These combinations
produce explicit errors.
The stored coefficient matrix \code{C} describes the latent equations;
\code{coefficientMap} maps latent coefficients to the original effect space.
Use \code{predict()}, \code{postPEV()}, and the fitted \code{uPevList} for
prediction uncertainty rather than directly indexing \code{C} by original
effect indices.

Earlier versions used \eqn{\Lambda\Lambda^{\mathsf T}+I_q}, estimating a
nonzero common specific variance on the fitted scale. To fit that former
covariance family, use \code{fam(x, k, variance="homogeneous")}.
}
\value{
A list containing the incidence/design matrix in \code{Z} and a compiled CovarianceFactor v2 descriptor in \code{covFactor}, for use inside \code{vsm}.
}
\examples{
\dontrun{
vsm(rrm(environment, k=2), ism(genotype))
# Former rrm() model: an estimated common specific variance
vsm(fam(environment, k=2, variance="homogeneous"), ism(genotype))
}
}
\seealso{
\code{\link{fam}}, \code{\link{ownm}}, \code{\link{vsm}}.
}

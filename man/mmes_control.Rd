\name{mmes_control}
\alias{mmes_control}
\title{Mixed Model Equation Solver Control and Experimental Options}
\description{
Constructs a structured configuration object for \code{\link{mmes}}, encapsulating
solver algorithms, memory budgets, iterative PCG parameters, and experimental numerical
paths without modifying global options.
}
\usage{
mmes_control(solver = c("auto", "ldlt", "cholmod", "pcg"),
             denseMemoryMB = NULL,
             experimental = list(),
             pcg = list(),
             acceleration = c("none", "aitken", "hybrid-ai"),
             factorScoreAugmentation = c("none", "fixed-shape", "profile"),
             recordDiagnostics = TRUE,
             verboseDiagnostics = FALSE,
             memoryPolicy = NULL,
             gpu = NULL,
             gpuMinSize = NULL,
             \dots)
}
\arguments{
  \item{solver}{Backend linear solver for coefficient matrix factorization.
  Options are \code{"auto"} (default; automatically selects \code{"cholmod"} when
  dense relationship precision is detected, otherwise \code{"ldlt"}),
  \code{"ldlt"} (simplicial LDLT), \code{"cholmod"} (supernodal blocked Cholesky),
  or \code{"pcg"} (preconditioned conjugate gradient with stochastic log-determinant).}

  \item{denseMemoryMB}{Upper bound on dense allocations (in megabytes) used by
  Schur-complement engines, block-chain/banded paths, and selected-inverse caches.
  Default is 1000 MB.}

  \item{memoryPolicy}{\code{"speed"} (default) stops with the memory needed when a fast path
  does not fit \code{denseMemoryMB}, instead of falling back to a slower path; \code{"limit"}
  keeps the fit within \code{denseMemoryMB} by accepting slower fallbacks. Currently enforced for
  the CHOLMOD selected-inverse caches, whose fallback (per-request solves) was measured several
  times slower, and opt-in exact LDLT trace panels. This is not a process-memory cap.}

  \item{experimental}{A named list of experimental solver paths. Each toggle accepts
  \code{"auto"} (policy-driven crossover), \code{"off"} (disabled), or \code{"force"}
  (enabled when structurally eligible):
  \describe{
    \item{\code{blockBanded}}{Exact block-banded elimination for assembled coefficient graphs with bandwidth \eqn{\le 4}.}
    \item{\code{factorSchurMultiFactor}}{Multi-factor latent Woodbury Schur complement engine.}
    \item{\code{dFirstML}}{D-first Schur complement elimination for maximum likelihood (ML) estimation.}
    \item{\code{matrixFreePcg}}{Matrix-free tensor operator for PCG random effect priors.}
    \item{\code{cholmodFinalSelected}}{Supernodal CHOLMOD Takahashi selected inverse for \code{computeCi=1}.}
    \item{\code{pdCertificates}}{Fast-path sufficient positive-definiteness certificates.}
    \item{\code{panelCi2}}{Batch panel size for full inverse computation in \code{computeCi=2} (default 64).}
    \item{\code{profile}}{\code{"force"} records per-phase timings and workload counters in \code{engineDiagnostics$profile}; default \code{"off"}.}
    \item{\code{blockSchur}}{Block elimination engines (dense Schur, chain, banded, latent): \code{"off"} uses the generic sparse factorization, \code{"force"} skips the density and group-size admission rules (memory limits still apply).}
    \item{\code{traceBatch}}{\code{"force"} batches exact LDLT random-score fallback solves in independent parallel panels; default \code{"off"}. No derivative entries are thresholded.}
    \item{\code{sparseFactors}}{\code{"force"} enables exact sparse AR(p)/ANTE precision operands, innovation-row ANTE covariance evaluation, and sparse SAR/CAR evaluation; default \code{"off"}. Dense covariance and derivative compatibility outputs remain available.}
    \item{\code{saddlePoint}}{\code{"force"} selects an experimental exact augmented MME
    factorization using pivoted sparse LU. Default and \code{"auto"} remain inactive.
    Gaussian Henderson REML/ML fits support diagonal, repeated-block, and explicit sparse
    residual covariance, including weights. Compact grid/ragged residuals, latent reduced rank,
    PQL, solve-only, rotation, acceleration, and factor-score augmentation are rejected.
    The deterministic solver setting still selects residual and final selected-inverse support;
    PCG is incompatible. Fitting avoids normal-equation products and retains exact traces and AI.
    Final \code{C} is assembled once for compatibility; this is not a C-free memory path.
    Mode 1 retains the Eigen selected-inverse subset, which can differ from the CHOLMOD subset.
    Diagnostics are in \code{engineDiagnostics$saddlePoint}. LU pivot fill is not bounded by
    \code{denseMemoryMB}; no speedup is guaranteed.}
  }}

  \item{pcg}{A named list of settings for the iterative PCG solver:
  \describe{
    \item{\code{tol}}{Convergence tolerance for iterative solves (default 1e-8).}
    \item{\code{maxIters}}{Maximum iterations per solve (0 = automatic, default 0).}
    \item{\code{traceProbes}}{Number of stochastic probes for Hutchinson traces (default 8).}
    \item{\code{lanczosSteps}}{Lanczos quadrature steps for stochastic log-determinant (default 20).}
    \item{\code{preconditioner}}{Preconditioner type: \code{"diagonal"},
      \code{"nystrom"}, or \code{"block-jacobi"}. These approximate the solve
      preconditioner, not the requested covariance model.}
    \item{\code{nystromRank}}{Rank for Nystrom preconditioner (default 32).}
    \item{\code{warmStart}}{Reuse compatible initial guesses (default \code{FALSE}).
      Every accepted solution must still satisfy the requested true-residual tolerance.}
    \item{\code{recycleRank}}{Maximum retained coarse-space dimension (default 0,
      disabled). Coarse products must be recomputed when the coefficient operator changes.}
    \item{\code{rhsMethod}}{\code{"independent"} (default) or opt-in \code{"block"}
      conjugate gradients for multiple right-hand sides.}
    \item{\code{blockSize}}{Maximum block-CG panel width (default 16).}
    \item{\code{residualRefresh}}{Periodic true-residual replacement interval
      (default 0, no periodic replacement). Final true-residual verification remains required.}
    \item{\code{logDetScaling}}{\code{"none"} (default) or \code{"diagonal"}.
      Diagonal scaling applies symmetric congruence and includes its exact determinant
      correction. Finite-probe estimates need not equal estimates without scaling.}
    \item{\code{traceEstimator}}{\code{"hutchinson"} (default) or experimental
      \code{"control-variate"} estimation. Estimator changes can change fitted estimates.}
    \item{\code{traceRank}}{Low-rank trace control-variate dimension (default 16).}
    \item{\code{accuracy}}{\code{"fixed"} (default) or \code{"adaptive"}.
      Adaptive mode doubles probe/depth budgets at the initial covariance point
      until empirical targets pass or their maxima are reached. It then freezes
      those budgets before accepting covariance updates, preserving a common
      stochastic objective throughout line search.}
    \item{\code{maxTraceProbes}, \code{maxLanczosSteps}}{Upper accuracy budgets.
      Defaults equal \code{traceProbes} and \code{lanczosSteps}; neither can be
      smaller than its initial setting.}
    \item{\code{logDetTol}, \code{scoreTol}}{Positive empirical accuracy targets
      (defaults 1e-3). These are not rigorous bounds on total likelihood or score error.}
    \item{\code{seed}}{Non-negative integer probe-stream identifier (default 0).
      The default preserves the existing deterministic probe stream.}
  }}

  \item{acceleration}{Convergence accelerator: \code{"none"} (default),
  \code{"aitken"}, or opt-in \code{"hybrid-ai"}. Hybrid AI uses Powell-damped
  BFGS curvature between exact average-information (AI) refreshes, at most five
  accepted iterations apart. Scores and likelihood acceptance remain unchanged.
  Between refreshes, the diagonal stabilizer follows the BFGS metric rather than
  exact AI. Fixed-parameter sensitivity solves are omitted on intermediate exact
  refreshes; full exact AI is restored before returning.

  Hybrid AI requires Gaussian identity-link Henderson fitting with \code{"ldlt"},
  \code{"cholmod"}, or \code{"auto"}. PCG, direct inversion, latent reduced-rank
  models, PQL, \code{solveOnly}, and factor-score augmentation are not supported.
  Fits with \code{vcc} ties retain exact AI and report \code{"constraint-map"} in
  \code{engineDiagnostics$hybridAI$fallbackReason}.

  Hybrid convergence requires an exact final relative proposal norm below
  \code{tolParConvNorm} and a projected, AI-diagonal-scaled score maximum no larger
  than \code{sqrt(max(tolParConvLL, tolParConvNorm))}. Reaching \code{nIters} or
  stalling without certification does not imply convergence. Final \code{InfMat}
  remains the exact AI/diagonal blend; \code{theta_se} always uses clean exact AI,
  independently of \code{computeCi}. Diagnostics report refreshes, BFGS updates,
  resets, certifications, actual sensitivity RHS columns, and skipped columns.
  This changes the optimization path, not the fitted covariance model; speedups
  depend on sensitivity costs and total iterations. Tolerances are never
  relaxed automatically; very tight requested step tolerances can return
  non-convergence when likelihood roundoff prevents further improvement.}

  \item{factorScoreAugmentation}{Factor-score augmentation: \code{"none"}, \code{"fixed-shape"}, or \code{"profile"}.}

  \item{recordDiagnostics}{Logical; whether to capture internal solver diagnostics in \code{fit$engineDiagnostics}. Default is \code{TRUE}.}

  \item{verboseDiagnostics}{Logical; whether to print diagnostic details during execution. Default is \code{FALSE}.}

  \item{gpu}{Logical; enable eligible CUDA dense kernels in the unified Henderson
  solver. Defaults to \code{getOption("sommer.mme.gpu", FALSE)}. Requires accessible
  NVIDIA hardware, compatible CUDA runtime, cuBLAS and cuSOLVER libraries, and a
  CUDA-enabled sommer build. Explicit requests stop when unavailable; CPU remains
  the default. Direct inversion and \code{solveOnly} remain CPU-only.}

  \item{gpuMinSize}{Positive integer; minimum matrix size for CUDA kernels.
  Defaults to \code{getOption("sommer.mme.gpuMinSize", 256L)}. Small eligible
  operations stay on CPU. GPU call counters appear in \code{fit$gpuDiagnostics}.}

  \item{\dots}{Additional named parameters or legacy experimental arguments.}
}
\details{
PCG-specific numerical options do not change the LDLT or CHOLMOD numerical
paths. They are opt-in; existing deterministic solver selection remains unchanged.
Unknown PCG names and invalid, non-scalar, non-finite or fractional integer
settings are errors, including settings merged into an existing engine object.

The PCG solve tolerance controls the relative linear-system residual. It does
not bound stochastic log-determinant, score, covariance-parameter, or uncertainty
errors. Probe standard errors describe Monte Carlo variability only, not
Lanczos truncation or total estimator bias. Compare final results with increased
probe counts and quadrature depths, and with direct solvers on manageable models.

REML PCG supports diagonal, model-aware block-Jacobi and operator-based Nystrom
preconditioners on assembled and eligible matrix-free coefficient operators.
Block-Jacobi groups random coefficients by relationship level and covariance
coordinate, with blocks capped by \code{blockSize}. Fixed coefficients use
bounded contiguous blocks. Every block must be positive definite.
\code{rhsMethod="independent"} retains packed scalar CG recurrences;
\code{"block"} uses coupled block CG with dependent-column deflation.
Numerical block breakdowns have diagnosed scalar fallbacks; invalid operators
and preconditioners remain errors. Aggressive periodic residual replacement can
increase iteration counts. Dense workspace budgets also govern PCG panels,
preconditioners, sketches and Lanczos storage.

Warm starts and coarse recycling retain only accepted-point coefficient and
probe solutions. Rejected line-search trials cannot replace accepted caches.
Every new covariance operator recomputes the recycled coarse Gram matrix, and
all returned solutions pass true-residual checks. Diagonally scaled SLQ does
not reuse untransformed Lanczos inverse guesses.

Trace control variates use an independently seeded low-rank space and an exact
projected correction. The same corrected estimator feeds all random and residual
score contractions, including weighted residual models. This identity does not
require positive-semidefinite covariance derivatives and does not guarantee
variance reduction for every model.

\code{fit$engineDiagnostics$pcg} reports operator calls and vector-equivalent
work, residual replacements, block fallbacks, warm/coarse reuse, probe estimates,
Monte Carlo standard errors and actual accuracy budgets. Adaptive final checks
require log-determinant probe standard error, nested half-depth sensitivity and
maximum free-score standard error below the requested targets. Failure emits a
warning and forces \code{convergence=FALSE}. Passing checks is empirical evidence,
not a rigorous error bound, exact REML certification or a guaranteed speedup.
Increasing budgets after initial calibration requires a new fit; budgets are not
silently changed within likelihood comparisons. ML and selected-inverse PCG
restrictions remain unchanged.

The known-covariance \code{solveOnly=TRUE} path uses its existing block-Jacobi
preconditioner. It additionally supports \code{warmStart} and
\code{residualRefresh}, but not REML estimator, recycling or block-RHS options.
Warm starts require a fitted \code{covPar} object with identical coefficient names;
the current operator and response always determine the verified final solution.
}
\value{
An object of class \code{"sommer_mmes_control"}.
}
\seealso{
\code{\link{mmes}}, \code{\link{generate_mme_workload}}, \code{\link{run_mme_benchmark}}
}
\examples{
# Default engine configuration
ctrl <- mmes_control()
print(ctrl)

# Force block-banded experimental path with custom memory budget
ctrl_banded <- mmes_control(
  solver = "cholmod",
  denseMemoryMB = 2000,
  experimental = list(blockBanded = "force")
)
print(ctrl_banded)
}

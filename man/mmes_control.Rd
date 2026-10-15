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
             acceleration = c("none", "aitken"),
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
  }}

  \item{pcg}{A named list of settings for the iterative PCG solver:
  \describe{
    \item{\code{tol}}{Convergence tolerance for iterative solves (default 1e-8).}
    \item{\code{maxIters}}{Maximum iterations per solve (0 = automatic, default 0).}
    \item{\code{traceProbes}}{Number of stochastic probes for Hutchinson traces (default 8).}
    \item{\code{lanczosSteps}}{Lanczos quadrature steps for stochastic log-determinant (default 20).}
    \item{\code{preconditioner}}{Preconditioner type: \code{"diagonal"} or \code{"nystrom"}.}
    \item{\code{nystromRank}}{Rank for Nystrom preconditioner (default 32).}
  }}

  \item{acceleration}{Convergence accelerator: \code{"none"} (default) or \code{"aitken"}.}

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

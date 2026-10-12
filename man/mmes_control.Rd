\name{mmes_control}
\alias{mmes_control}
\alias{mme_control}
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
             \dots)

mme_control(\dots)
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

  \item{\dots}{Additional named parameters or legacy experimental arguments.}
}
\value{
An object of class \code{"sommer_mme_control"}.
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

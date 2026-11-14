\name{mmes_control}
\alias{mmes_control}
\title{

Numerical controls for mmes

}
\description{

Creates the validated object passed to the \code{control} argument of
\code{\link{mmes}}. It collects every \emph{numerical} choice of a fit:
the equation formulation, the linear solvers, covariance updates,
convergence, output, memory/threads and diagnostics. Statistical choices
(REML vs ML, exact vs Laplace vs PQL integration, fit vs solve) do not belong
here; they are declared in \code{\link{mmes_method}}. Starting values for
nonlinear parameters are given in \code{mmes(start=)} and family nuisance
parameters in \code{mmes(familyParameters=)}.

Each argument is a module and is given as a named list; only the fields you
supply are changed and all others keep their defaults. For example,
\code{mmes_control(solve=list(method="pcg", tol=1e-6), execution=list(seed=1))}.
Unknown modules or fields, invalid values, and combinations that are not
mathematically admissible (for example PCG fields without
\code{solve$method="pcg"}) stop with an informative error that names the
offending fields. A request can also be valid on its own but ineligible for
the model being fitted; \code{mmes} then stops rather than silently changing
the algorithm. The defaults reproduce the established sommer algorithm, so
\code{mmes_control()} is the right choice for most analyses.

}
\usage{

mmes_control(formulation = "henderson", optimizer = list(), mode = list(),
    working = list(), likelihood = list(), population = list(),
    inference = list(), representation = list(), solve = list(),
    logdet = list(), score = list(), curvature = list(), update = list(),
    proposal = list(), convergence = list(), output = list(),
    execution = list(), diagnostics = list())

}
\arguments{

  \item{formulation}{A single string selecting the equation domain.
  \itemize{
    \item \code{"henderson"} (default): coefficient-space mixed model
    equations (fixed and random effects together). Sparse and fast when the
    number of coefficients is moderate or the design is sparse. Required by
    most advanced options (PCG, Schur/banded elimination, compression,
    factor scores, selected inverse, PQL and Laplace GLMMs).
    \item \code{"direct"}: observation-space direct inversion of the
    marginal covariance. Can be preferable when there are many more random
    coefficients than observations. Only models that support it accept it.
  }
  The factorization backend is chosen separately in \code{solve$method}.}

  \item{optimizer}{Outer optimizer for nonlinear-parameter and native Laplace
  likelihood fits (Gaussian covariance engines use the AI/update modules
  below instead). Fields:
  \itemize{
    \item \code{maxit}: \code{NULL} (default, engine default) or a positive
    integer limiting outer iterations. Do not also set a different
    \code{convergence$maxIters}; the limit must be specified once.
    \item \code{gradient}: \code{"auto"} (default) uses native analytic
    gradients where available and reported numerical derivatives for
    reference expressions; \code{"analytic"} requires analytic scores and
    errors if unavailable; \code{"central"} uses central differences.
    \item \code{algorithm}: only \code{"nlminb"} is currently available.
    \item \code{scoreTol}: \code{NULL} (default, 0.001) or a positive number,
    the outer-score certificate for nonlinear specification fits. Not allowed
    for non-nonlinear models.
  }}

  \item{mode}{Conditional random-effect mode search used by Laplace and
  nonlinear fits (not by Gaussian or PQL working fits). Fields:
  \itemize{
    \item \code{maxit}: maximum Newton iterations per mode (default
    \code{100}).
    \item \code{tol}: score tolerance for mode convergence (default
    \code{1e-8}).
    \item \code{stepTol}: step tolerance that detects stalled refinement
    (default \code{1e-8}); a non-default value is supported only by
    nonlinear modes, native GLMM modes use \code{tol}.
    \item \code{warmStart}: \code{"auto"} (default; reuse successful modes
    for native GLMMs, deterministic cold starts for nonlinear models),
    \code{"off"} or \code{"on"} (errors for nonlinear models).
    \item \code{starts}: \code{"deterministic"} (only choice).
  }}

  \item{working}{Outer iterations of the penalized quasi-likelihood (PQL)
  working fit used for non-Gaussian mixed models by default. Fields:
  \code{maxit} (default \code{20}), \code{tol} (default \code{1e-5}),
  \code{estimateTheta} (default \code{FALSE}; \code{TRUE} estimates the
  negative binomial theta conditionally on the random effects and requires
  \code{MASS::negative.binomial}), and \code{warmStart} (default \code{TRUE},
  start each outer iteration from the previous solution). The inner Gaussian
  fits are controlled by \code{convergence}.}

  \item{likelihood}{Algebra used by the native Laplace likelihood (see
  \code{mmes_method(integration="laplace")}); it does not choose the Gaussian
  solver. Fields:
  \itemize{
    \item \code{backend}: \code{"auto"} (default), \code{"random"} (algebra
    in random-effect space) or \code{"observation"} (in observation space).
    All compute the same approximation; only \code{"auto"} is accepted for
    non-Laplace fits.
    \item \code{trace}: \code{"panels"} (default; bounded dense solve panels
    sized by \code{execution$denseMemoryMB}) or \code{"selected"}.
    \item \code{blocks}: \code{"off"} (default) or \code{"on"} to use
    block-structured post-fit operations.
  }}

  \item{population}{Gauss quadrature for response-scale population
  predictions of nonlinear models. Fields: \code{order} (initial order,
  default \code{20}), \code{maxOrder} (maximum order, default \code{160},
  at most 512 and at least twice \code{order}) and \code{tol} (accuracy,
  default \code{1e-7}). If the accuracy is not met at \code{maxOrder} or the
  workspace limit is hit, an error is raised.}

  \item{inference}{Post-fit numerical uncertainty for nonlinear models.
  Fields: \code{step} (relative central-difference step, default
  \code{1e-4}, must be below 0.1), \code{derivativeTol} (derivative
  acceptance tolerance, default \code{5e-3}, below 1), \code{curvatureTol}
  (curvature tolerance, default \code{1e-8}, below 1) and
  \code{maxCondition} (largest condition number accepted for the Hessian,
  default \code{1e10}, above 1). Standard errors are reported as \code{NA}
  when these safeguards fail.}

  \item{representation}{How equations and covariance operators are
  represented. Fields:
  \itemize{
    \item \code{compression}: \code{"off"} (default), \code{"groups"}
    (collapse exactly identical observation groups) or \code{"classes"}
    (also cache cross-products and derivative operands). Exact reductions
    for Gaussian Henderson fits with independent residuals and fixed
    diagonal weights.
    \item \code{factorScores}: \code{"none"} (default), \code{"fixed-shape"}
    or \code{"profile"} for equivalent latent representations of eligible
    low-rank-plus-diagonal (factor analytic/reduced rank) models; shape
    standard errors are unavailable with \code{"profile"}.
    \item \code{fixedRank}: \code{"error"} (default; stop if the fixed design
    is rank deficient) or \code{"drop"} (fit a full-rank subset using sparse
    QR; omitted coefficients are returned as zeros, flagged in
    \code{summary(fit)$betas}, and are identification constraints, not
    estimates).
    \item \code{equations}: \code{"normal"} (default), \code{"saddle-point"}
    (augmented system; avoids forming the compatibility coefficient matrix
    when \code{output$coefficientMatrix=FALSE}) or the experimental
    \code{"conditional"} (reusable sparse latent equations).
    \item \code{sparseFactors}: \code{"auto"}, \code{"off"} or \code{"force"}
    (\code{TRUE}/\code{FALSE} accepted). Uses exact sparse precision
    operands for autoregressive, antedependence, SAR and CAR factors. Auto
    also keeps complete-grid separable residuals banded.
    \item \code{gridAugment}: \code{"auto"}, \code{"off"} or \code{"force"}.
    Adds one zero-response row and dummy fixed effect per missing cell of a
    separable gridded residual so blocks are complete banded Kronecker
    precisions. Likelihood, variances and BLUEs/BLUPs are unchanged; outputs
    are mapped back to the observed design.
    \item \code{matrixFree}: \code{"auto"} (default), \code{"off"} or
    \code{"force"}; avoids assembling full Kronecker priors during PCG.
  }}

  \item{solve}{Linear solves and exact eliminations. Fields:
  \itemize{
    \item \code{method}: \code{"auto"} (default; the model decides),
    \code{"ldlt"} (sparse LDLT), \code{"cholmod"} (supernodal sparse Cholesky)
    or \code{"pcg"} (preconditioned conjugate gradients with stochastic
    traces; for very large sparse systems).
    \item PCG fields (only valid with \code{method="pcg"} or known-covariance
    solving where noted): \code{tol} (relative tolerance, default
    \code{1e-8}); \code{maxIters} (CG limit, \code{0} = automatic);
    \code{preconditioner} (\code{"diagonal"} default, \code{"nystrom"},
    \code{"block-jacobi"}; solving permits diagonal and block-jacobi only);
    \code{nystromRank} (default \code{32}); \code{warmStart} (reuse accepted
    solutions, default \code{FALSE}); \code{recycleRank} (accepted points
    reused, \code{0} disables, default \code{0}); \code{recycleMethod}
    (\code{"coarse"} default or \code{"krylov"}, which needs positive
    \code{recycleRank}); \code{recycleState} and \code{exportRecycle} (import
    or export a serializable Krylov basis across compatible fits; default
    \code{NULL}/\code{FALSE}); \code{rhsMethod} (\code{"independent"} default
    or \code{"block"}); \code{blockSize} (block width, default \code{16});
    \code{residualRefresh} (periodic residual replacement, \code{0}
    disables).
    \item Exact elimination switches, each \code{"auto"} (default),
    \code{"off"} or \code{"force"}: \code{blockSchur} (dense Schur, chain,
    banded or latent elimination), \code{blockBanded} (block-banded
    elimination, bandwidth at most 4), \code{factorSchurMultiFactor}
    (multi-factor Woodbury/Schur for latent low-rank layouts) and
    \code{dFirstML} (eliminate the random block before the fixed-effect Schur
    complement in ML). A forced request still obeys structural and memory
    checks.
  }}

  \item{logdet}{Log-determinant approximation used by PCG. Fields:
  \code{scaling} (\code{"none"} default or \code{"diagonal"}, a symmetric
  scaling with exact determinant correction) and \code{lanczosSteps}
  (stochastic Lanczos quadrature depth, integer \eqn{\ge 2}, default
  \code{20}).}

  \item{score}{Covariance-parameter derivatives and trace estimation.
  Fields: \code{traceProbes} (probe budget for stochastic traces, default
  \code{8}), \code{estimator} (\code{"hutchinson"} default or
  \code{"control-variate"}, which uses a low-rank correction of dimension
  \code{rank}, default \code{16}); the previous three apply only to PCG.
  Also \code{covarianceChecks} (\code{"auto"} default, \code{"off"},
  \code{"force"}; skips redundant positive-definiteness checks when
  certified), \code{traceBatch} (\code{"off"} default, \code{"auto"},
  \code{"force"}; batch exact LDLT traces) and \code{denseCrossproduct}
  (\code{"auto"} default; uses an exact weighted sparse Gram product for
  dense designs).}

  \item{curvature}{Information matrix used by the optimizer. \code{method}:
  \code{"exact-ai"} (default; exact average information) or
  \code{"hybrid-ai"} (damped BFGS between exact refreshes).
  \code{assembly}: \code{"sensitivity"} (default) or
  \code{"augmented-working"} (exact coefficient-space working-variate
  absorption for Gaussian Henderson fits). \code{reduced}: \code{"off"}
  (default), \code{"auto"} or \code{"force"}; exact independent-direction
  AI assembly when \code{vcc} reduces the parameters.}

  \item{update}{How covariance parameters are updated. \code{method}:
  \code{"ai"} (default average-information Newton steps), \code{"em"}
  (conditional EM for eligible factor analytic models) or
  \code{"conditional-fa"} (legacy EM selector). \code{recovery}: \code{"em"}
  (default) or \code{"pxem"}, the EM-type direction used for recovery.
  \code{blendWeight}: \code{NULL} (default), a non-negative scalar or vector
  weighting the recovery direction (an integer larger than 1 requests a
  decaying schedule). \code{commonScale}: \code{"off"} (default),
  \code{"auto"} or \code{"force"}; analytic common variance scale plus
  profiled-shape AI for Gaussian fits where all product scales are free.}

  \item{proposal}{Scaling of proposed steps. \code{stepWeight}: \code{NULL}
  (default) or positive scalar/vector multiplying the step; likelihood
  safeguards stay active. \code{extrapolation}: \code{"none"} (default) or
  \code{"aitken"} (safeguarded Aitken extrapolation).}

  \item{convergence}{Inner iterations and stopping rules. Fields:
  \code{maxIters} (maximum covariance iterations, default \code{30}; reaching
  it is not proof of convergence), \code{likelihoodTol} and
  \code{parameterTol} (default \code{1e-4}; relative change in log-likelihood
  and in parameters), and, for PCG: \code{accuracy} (\code{"fixed"} default
  or \code{"adaptive"}, which grows probe budgets up to \code{maxTraceProbes}
  and \code{maxLanczosSteps}, defaulting to \code{score$traceProbes} and
  \code{logdet$lanczosSteps} and never smaller than them), \code{logDetTol}
  and \code{scoreTol} (empirical accuracy targets, default \code{0.001}).}

  \item{output}{What is computed after the fit. Fields:
  \itemize{
    \item \code{inverse}: \code{0} (default; no coefficient inverse or PEVs),
    \code{1} (random-effect PEVs from the selected inverse, without full
    \code{Ci}) or \code{2} (full \code{Ci} and PEVs; only advisable for small
    systems). Defer this work with \code{0} and call
    \code{\link{postPEV}} later. Native Laplace fits require \code{0}.
    \item \code{covariance}: \code{"full"} (default; standard errors of
    covariance parameters) or \code{"none"} (skip final covariance
    uncertainty; point estimates remain available).
    \item \code{coefficientMatrix}: \code{TRUE} (default) or \code{FALSE}
    (skip assembling the compatibility matrix \code{C}; needs saddle-point
    equations and \code{inverse=0}; \code{C}-based methods then need a refit).
    \item \code{setup}: \code{FALSE} (default) or \code{TRUE} to return the
    prepared design, covariance descriptors and \code{vcParams} table without
    fitting (useful for building \code{vcc} constraints).
    \item \code{selectedInverse}: \code{"auto"} (default), \code{"off"} or
    \code{"force"}; exact selected inverse from CHOLMOD with
    \code{inverse=1}.
    \item \code{panelSize}: width of inverse-solve batches (default
    \code{64}), also used by \code{postPEV} and Laplace predictions.
  }}

  \item{execution}{Resources and reproducibility. Fields:
  \code{denseMemoryMB} (admission budget for dense workspaces, default
  \code{1000}; not a process memory cap), \code{memoryPolicy}
  (\code{"speed"} default or \code{"limit"} to enforce the budget),
  \code{gpu} (\code{FALSE} default; \code{TRUE} uses eligible dense kernels on
  CUDA builds with usable NVIDIA libraries), \code{gpuMinSize} (minimum
  kernel size for the GPU, default \code{256}), \code{blasThreads}
  (\code{NULL} default, a positive integer or \code{"auto"}; BLAS only,
  OpenMP is configured separately), \code{inverseTol} (numerical
  stabilization and variance floor, default \code{1e-6}; not the PCG
  tolerance) and \code{seed} (non-negative integer identifying the stochastic
  probe stream, default \code{0}).}

  \item{diagnostics}{Messages and profiling. Fields: \code{verbose}
  (\code{TRUE} default; print the iteration table and fit messages),
  \code{dateWarning} (\code{TRUE} default; package-age warning) and
  \code{profile} (\code{"off"} default, \code{"auto"} or \code{"force"};
  phase timings stored in \code{engineDiagnostics$profile}, and see
  \code{\link{profile_nonlinear}} for nonlinear fits).}

}
\details{

\strong{Switch-like fields.} Fields documented as \code{"auto"}, \code{"off"},
\code{"force"} also accept \code{TRUE} (force) and \code{FALSE} (off).
\code{"auto"} lets sommer use the option when the model is eligible, while
\code{"force"} requests it and errors if it is not admissible. A forced
request never bypasses mathematical, structural or memory checks.

\strong{Global options.} Experimental-mode defaults, the memory budget and
GPU settings can be set through R options, e.g.
\code{options(sommer.mme.denseMemoryMB=2000)}. Values given explicitly in
\code{mmes_control()} always override option-derived defaults, so use them
for reproducible analyses.

\strong{Common recipes.}
\itemize{
  \item Quiet fit: \code{mmes_control(diagnostics=list(verbose=FALSE))}.
  \item More iterations: \code{mmes_control(convergence=list(maxIters=100))}.
  \item Random-effect PEVs at fit time: \code{mmes_control(output=list(inverse=1))};
  or leave \code{inverse=0} and call \code{\link{postPEV}}.
  \item Skip standard errors of variance components for speed:
  \code{mmes_control(output=list(covariance="none"))}.
  \item Very large sparse system:
  \code{mmes_control(solve=list(method="pcg"), execution=list(seed=1))}.
  \item Rank-deficient fixed design:
  \code{mmes_control(representation=list(fixedRank="drop"))}.
  \item Inspect the model without fitting:
  \code{mmes_control(output=list(setup=TRUE))}.
  \item Limit threads and memory:
  \code{mmes_control(execution=list(blasThreads=2, denseMemoryMB=500))}.
}

\strong{Migration.} Numerical arguments removed from \code{mmes()} are now
fields of the modules above, e.g. \code{henderson} is
\code{formulation}, \code{solver} is \code{solve$method}, \code{nIters} is
\code{convergence$maxIters}, \code{tolParConvLL} and \code{tolParConvNorm}
are \code{convergence$likelihoodTol} and \code{convergence$parameterTol},
\code{tolParInv} is \code{execution$inverseTol}, \code{computeCi} is
\code{output$inverse}, \code{returnParam} is \code{output$setup},
\code{pqlControl} is \code{working}, and \code{verbose}/\code{dateWarning}
are \code{diagnostics} fields. \code{REML} and \code{solveOnly} moved to
\code{\link{mmes_method}}. The complete table is in the
\code{sommer.mmes.options} vignette.

}
\value{

A list of class \code{sommer_mmes_control} with one element per module.
It is validated again by \code{mmes()} and stored in the fitted object as
\code{fit$control}. \code{mmes()} accepts only objects created by this
function.

}
\references{

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits using
the R package sommer. PLoS ONE 2016, 11(6): doi:10.1371/journal.pone.0156744

Gilmour et al. 1995. Average Information REML: An efficient algorithm for
variance parameter estimation in linear mixed models. Biometrics 51(4):1440-1450.

}
\seealso{

\code{\link{mmes}}, \code{\link{mmes_method}}, \code{\link{postPEV}},
\code{\link{vcc}}, \code{\link{profile_nonlinear}}; the vignette
\code{vignette("sommer.mmes.options")} for the full field reference, defaults,
eligibility rules and decision graph, and \code{vignette("sommer.machinery")}
for the engine.

}
\examples{

dat <- data.frame(x=seq_len(12), y=1+seq_len(12)/2+sin(seq_len(12))/10)

## defaults, silenced
ctl <- mmes_control(diagnostics=list(verbose=FALSE, dateWarning=FALSE))
fit <- mmes(y~x, data=dat, control=ctl)
coef(fit)

## modules only override the fields supplied
ctl2 <- mmes_control(convergence=list(maxIters=50),
                     output=list(covariance="none"),
                     diagnostics=list(verbose=FALSE, dateWarning=FALSE))
ctl2$convergence$maxIters
ctl2$convergence$likelihoodTol

## invalid choices are rejected
try(mmes_control(solve=list(method="gauss")))

}

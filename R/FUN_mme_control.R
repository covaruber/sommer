# Structured solver engine configuration for mmes()
# Encapsulates solver selection, memory budgets, iterative (PCG) parameters,
# and experimental path toggles without polluting global options.

.mme_experimental_mode <- function(val, default="auto", optName=NULL){
  if(is.null(val)){
    if(!is.null(optName)){
      optVal <- getOption(optName, NULL)
      if(!is.null(optVal)){
        if(is.logical(optVal) && length(optVal) == 1L && !is.na(optVal)){
          return(if(optVal) "force" else "off")
        }
        if(is.character(optVal) && length(optVal) == 1L){
          val <- optVal
        } else {
          stop(sprintf("%s must be TRUE or FALSE.", optName), call.=FALSE)
        }
      }
    }
    if(is.null(val)) return(default)
  }
  if(is.logical(val) && length(val) == 1L){
    if(is.na(val)){
      if(!is.null(optName)){
        stop(sprintf("%s must be TRUE or FALSE.", optName), call.=FALSE)
      }
      return(default)
    }
    return(if(val) "force" else "off")
  }
  if(is.character(val) && length(val) == 1L){
    valLower <- tolower(val)
    if(valLower %in% c("auto", "off", "force", "true", "false")){
      if(valLower == "true") return("force")
      if(valLower == "false") return("off")
      return(valLower)
    }
  }
  stop(sprintf("Experimental toggle must be 'auto', 'off', 'force', TRUE, or FALSE (got '%s').",
               as.character(val)), call.=FALSE)
}

.mme_pcg_control <- function(pcg = list(), checkLimits = TRUE){
  defaults <- list(
    tol = 1e-8, maxIters = 0L, traceProbes = 8L, lanczosSteps = 20L,
    preconditioner = "diagonal", nystromRank = 32L, warmStart = FALSE,
    recycleRank = 0L, rhsMethod = "independent", blockSize = 16L,
    residualRefresh = 0L, logDetScaling = "none",
    traceEstimator = "hutchinson", traceRank = 16L, accuracy = "fixed",
    maxTraceProbes = NULL, maxLanczosSteps = NULL,
    logDetTol = 1e-3, scoreTol = 1e-3, seed = 0L
  )
  if(!is.list(pcg) || (length(pcg) &&
     (is.null(names(pcg)) || anyNA(names(pcg)) || any(!nzchar(names(pcg))) ||
      anyDuplicated(names(pcg))))){
    stop("pcg must be a list with unique, non-empty names.", call.=FALSE)
  }
  unknown <- setdiff(names(pcg), names(defaults))
  if(length(unknown)){
    stop("Unknown pcg settings: ", paste(unknown, collapse=", "), ".", call.=FALSE)
  }
  out <- utils::modifyList(defaults, pcg, keep.null=TRUE)
  if(is.null(out$maxTraceProbes)) out$maxTraceProbes <- out$traceProbes
  if(is.null(out$maxLanczosSteps)) out$maxLanczosSteps <- out$lanczosSteps
  for(name in c("tol", "logDetTol", "scoreTol")){
    value <- out[[name]]
    if(!is.numeric(value) || is.complex(value) || length(value) != 1L ||
       !is.finite(value) || value <= 0){
      stop("pcg$", name, " must be one positive finite number.", call.=FALSE)
    }
  }
  minima <- c(maxIters=0L, traceProbes=1L, lanczosSteps=2L, nystromRank=1L,
              recycleRank=0L, blockSize=1L, residualRefresh=0L, traceRank=1L,
              maxTraceProbes=1L, maxLanczosSteps=2L, seed=0L)
  for(name in names(minima)){
    value <- out[[name]]
    if(!is.numeric(value) || is.complex(value) || length(value) != 1L ||
       !is.finite(value) || value < minima[[name]] ||
       value > .Machine$integer.max || value != floor(value)){
      stop("pcg$", name, " must be one integer between ", minima[[name]],
           " and ", .Machine$integer.max, ".", call.=FALSE)
    }
    out[[name]] <- as.integer(value)
  }
  choices <- list(
    preconditioner=c("diagonal", "nystrom", "block-jacobi"),
    rhsMethod=c("independent", "block"),
    logDetScaling=c("none", "diagonal"),
    traceEstimator=c("hutchinson", "control-variate"),
    accuracy=c("fixed", "adaptive")
  )
  for(name in names(choices)){
    value <- out[[name]]
    if(!is.character(value) || length(value) != 1L || is.na(value) ||
       !tolower(value) %in% choices[[name]]){
      stop("pcg$", name, " must be one of ",
           paste(shQuote(choices[[name]]), collapse=", "), ".", call.=FALSE)
    }
    out[[name]] <- tolower(value)
  }
  if(!is.logical(out$warmStart) || length(out$warmStart) != 1L ||
     is.na(out$warmStart)){
    stop("pcg$warmStart must be a single TRUE/FALSE value.", call.=FALSE)
  }
  if(checkLimits && (out$maxTraceProbes < out$traceProbes ||
     out$maxLanczosSteps < out$lanczosSteps)){
    stop("PCG accuracy limits must not be smaller than traceProbes and lanczosSteps.",
         call.=FALSE)
  }
  out
}

#' Mixed Model Equation Engine Control
#'
#' Builds a structured control object defining solver algorithms, memory budgets,
#' iterative settings, and experimental numerical paths for \code{\link{mmes}}.
#'
#' @param solver Backend solver for coefficient matrix factorization. One of \code{"auto"},
#'   \code{"ldlt"}, \code{"cholmod"}, or \code{"pcg"}.
#' @param denseMemoryMB Upper bound on dense memory allocations (in megabytes) used by
#'   Schur-complement engines, block-chain / banded paths, and selected-inverse caches. Default is 1000.
#' @param memoryPolicy \code{"speed"} (default) stops with the memory needed when a fast path
#'   does not fit \code{denseMemoryMB}, instead of falling back to a slower path; \code{"limit"}
#'   keeps the fit within \code{denseMemoryMB} by accepting slower fallbacks. Currently enforced for
#'   the CHOLMOD selected-inverse caches, whose fallback (per-request solves) was measured several
#'   times slower, and opt-in exact LDLT trace panels. This is not a process-memory cap.
#' @param experimental A named list of experimental path controls. Each toggle accepts \code{"auto"},
#'   \code{"off"}, or \code{"force"} (or boolean \code{TRUE}/\code{FALSE}):
#'   \itemize{
#'     \item \code{blockBanded}: Exact block-banded CHOLMOD elimination for banded graphs.
#'     \item \code{factorSchurMultiFactor}: Multi-factor latent Woodbury Schur complement engine.
#'     \item \code{dFirstML}: D-first Schur elimination for ML estimation.
#'     \item \code{matrixFreePcg}: Matrix-free tensor operator for PCG priors.
#'     \item \code{cholmodFinalSelected}: Supernodal CHOLMOD Takahashi selected inverse for mode 1.
#'     \item \code{pdCertificates}: Sufficient positive-definiteness certificates fast-path.
#'     \item \code{panelCi2}: Batch panel size for full inverse in \code{computeCi=2} (default 32).
#'     \item \code{profile}: \code{"force"} records per-phase timings and workload counters in
#'       \code{engineDiagnostics$profile}; default \code{"off"}.
#'     \item \code{blockSchur}: block elimination engines (dense Schur, chain, banded, latent).
#'       \code{"off"} uses the generic sparse factorization; \code{"force"} skips the density and
#'       group-size admission rules (memory limits still apply).
#'     \item \code{traceBatch}: \code{"force"} batches exact LDLT random-score fallback solves
#'       in independent parallel panels; default \code{"off"}. No derivative entries are thresholded.
#'     \item \code{sparseFactors}: \code{"force"} enables exact sparse AR(p)/ANTE precision operands,
#'       innovation-row ANTE covariance evaluation, and sparse SAR/CAR evaluation; default \code{"off"}.
#'     \item \code{saddlePoint}: \code{"force"} selects an experimental exact augmented MME
#'       factorization using pivoted sparse LU. Default and \code{"auto"} remain inactive.
#'       Gaussian Henderson REML/ML fits support diagonal, repeated-block, and explicit sparse
#'       residual covariance, including weights. Compact grid/ragged residuals, latent reduced
#'       rank, PQL, solve-only, rotation, acceleration, and factor-score augmentation are rejected.
#'       The deterministic solver setting still selects residual and final selected-inverse
#'       support; PCG is incompatible. Fitting avoids normal-equation products and uses exact
#'       traces and AI. Final \code{C} is assembled once for compatibility, so this is not a
#'       C-free memory path. Mode 1 retains the Eigen selected-inverse subset, which can differ
#'       from the CHOLMOD subset. Diagnostics are in \code{engineDiagnostics$saddlePoint}.
#'       LU pivot fill is not bounded by \code{denseMemoryMB}; no speedup is guaranteed.
#'   }
#' @param pcg A named list of iterative PCG solver settings:
#'   \itemize{
#'     \item \code{tol}: Convergence tolerance for iterative CG solves (default 1e-8).
#'     \item \code{maxIters}: Maximum iterations for PCG solves (0 = automatic, default 0).
#'     \item \code{traceProbes}: Number of stochastic probes for Hutchinson traces (default 8).
#'     \item \code{lanczosSteps}: Lanczos quadrature steps for stochastic log determinant (default 20).
#'     \item \code{preconditioner}: \code{"diagonal"}, \code{"nystrom"}, or
#'       \code{"block-jacobi"}.
#'     \item \code{nystromRank}: Target rank for Nystrom preconditioner (default 32).
#'     \item \code{warmStart}, \code{recycleRank}: Reuse compatible solve guesses
#'       and a bounded coarse space (defaults FALSE and 0).
#'     \item \code{rhsMethod}, \code{blockSize}: Independent or block CG and its
#'       panel width (defaults \code{"independent"} and 16).
#'     \item \code{residualRefresh}: Recursive-residual replacement interval;
#'       zero disables periodic replacement. Final true residuals are checked.
#'     \item \code{logDetScaling}: \code{"none"} or symmetric \code{"diagonal"}
#'       congruence for stochastic log determinants.
#'     \item \code{traceEstimator}, \code{traceRank}: Hutchinson or low-rank
#'       \code{"control-variate"} trace estimation and its rank (default 16).
#'     \item \code{accuracy}: \code{"fixed"} or opt-in \code{"adaptive"} accuracy.
#'       Adaptive mode doubles probe/depth budgets at the initial covariance point,
#'       then freezes the objective before line search. Final empirical checks
#'       can refuse convergence; they do not certify total statistical error.
#'     \item \code{maxTraceProbes}, \code{maxLanczosSteps}: Accuracy limits;
#'       default to the initial probe count and Lanczos depth.
#'     \item \code{logDetTol}, \code{scoreTol}: Positive empirical accuracy targets
#'       (defaults 1e-3); these are not rigorous total-error bounds.
#'     \item \code{seed}: Non-negative integer probe stream identifier (default 0).
#'   }
#' @param acceleration Convergence accelerator. One of \code{"none"}, \code{"aitken"},
#'   or \code{"hybrid-ai"}. The opt-in hybrid uses damped BFGS between exact AI
#'   refreshes and always returns exact final covariance information. It supports
#'   Gaussian Henderson fits with deterministic solvers; tied parameters use exact AI.
#'   Exact refreshes occur at most five accepted iterations apart and for final
#'   certification. Between refreshes, diagonal stabilization uses the current
#'   BFGS metric. PCG, direct inversion, latent RR, PQL, solve-only and factor-score
#'   augmentation are not supported. Diagnostics in \code{engineDiagnostics$hybridAI}
#'   report refreshes, resets, certifications and actual/skipped sensitivity columns.
#'   Hybrid convergence requires an exact relative proposal norm below
#'   \code{tolParConvNorm} and an AI-diagonal-scaled projected score maximum no larger
#'   than \code{sqrt(max(tolParConvLL, tolParConvNorm))}. Tolerances are never relaxed
#'   automatically. Final \code{InfMat} retains its exact AI/diagonal blend;
#'   \code{theta_se} uses clean exact AI independently of \code{computeCi}.
#' @param factorScoreAugmentation Factor-score augmentation strategy. One of \code{"none"},
#'   \code{"fixed-shape"}, or \code{"profile"}.
#' @param recordDiagnostics Logical; whether to capture detailed engine diagnostics in the fitted object. Default is TRUE.
#' @param verboseDiagnostics Logical; whether to print diagnostic details during execution. Default is FALSE.
#' @param gpu Logical; enable eligible CUDA dense kernels. Defaults to
#'   \code{getOption("sommer.mme.gpu", FALSE)}. Requires usable NVIDIA hardware,
#'   CUDA libraries, and a CUDA-enabled package build.
#' @param gpuMinSize Positive integer; minimum matrix size for CUDA kernels.
#'   Defaults to \code{getOption("sommer.mme.gpuMinSize", 256L)}.
#' @param ... Additional experimental toggles or legacy options passed as named arguments.
#'
#' @return An object of class \code{"sommer_mmes_control"}.
#' @export
mmes_control <- function(solver = c("auto", "ldlt", "cholmod", "pcg"),
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
                         ...){
  dots <- list(...)
  if(is.null(gpu)) gpu <- getOption("sommer.mme.gpu", FALSE)
  if(!is.logical(gpu) || length(gpu) != 1L || is.na(gpu)){
    stop("gpu must be a single TRUE/FALSE value.", call.=FALSE)
  }
  if(is.null(gpuMinSize)) gpuMinSize <- getOption("sommer.mme.gpuMinSize", 256L)
  if(!is.numeric(gpuMinSize) || length(gpuMinSize) != 1L ||
     !is.finite(gpuMinSize) || gpuMinSize < 1 || gpuMinSize > .Machine$integer.max ||
     gpuMinSize != floor(gpuMinSize)){
    stop("gpuMinSize must be a positive integer.", call.=FALSE)
  }
  if(missing(solver) || is.null(solver)){
    solver <- "auto"
  } else {
    solverChoices <- c("auto", "ldlt", "cholmod", "pcg")
    if(!is.character(solver) || length(solver) < 1L || is.na(solver[1L]) ||
       !tolower(solver[1L]) %in% solverChoices){
      stop("solver must be one of 'auto', 'ldlt', 'cholmod', or 'pcg'.", call.=FALSE)
    }
    solver <- match.arg(tolower(solver[1L]), c("auto", "ldlt", "cholmod", "pcg"))
  }

  if(missing(acceleration) || is.null(acceleration)){
    acceleration <- "none"
  } else {
    acceleration <- match.arg(tolower(acceleration[1L]), c("none", "aitken", "hybrid-ai"))
  }

  if(missing(factorScoreAugmentation) || is.null(factorScoreAugmentation)){
    factorScoreAugmentation <- "none"
  } else {
    factorScoreAugmentation <- match.arg(tolower(factorScoreAugmentation[1L]),
                                         c("none", "fixed-shape", "profile"))
  }

  # Resolve denseMemoryMB with fallback to getOption("sommer.mme.denseMemoryMB", 1000)
  if(is.null(denseMemoryMB)){
    optMB <- getOption("sommer.mme.denseMemoryMB", NULL)
    if(!is.null(optMB)){
      if(!is.numeric(optMB) || length(optMB) != 1L || !is.finite(optMB) || optMB <= 0){
        stop("sommer.mme.denseMemoryMB must be a positive finite number.", call.=FALSE)
      }
      denseMemoryMB <- as.double(optMB)
    } else {
      denseMemoryMB <- 1000.0
    }
  } else {
    if(!is.numeric(denseMemoryMB) || length(denseMemoryMB) != 1L ||
       !is.finite(denseMemoryMB) || denseMemoryMB <= 0){
      stop("sommer.mme.denseMemoryMB must be a positive finite number.", call.=FALSE)
    }
    denseMemoryMB <- as.double(denseMemoryMB)
  }

  if(is.null(memoryPolicy)) memoryPolicy <- getOption("sommer.mme.memoryPolicy", "speed")
  if(!is.character(memoryPolicy) || length(memoryPolicy) != 1L ||
     !tolower(memoryPolicy) %in% c("speed", "limit")){
    stop("memoryPolicy must be 'speed' or 'limit'.", call.=FALSE)
  }
  memoryPolicy <- tolower(memoryPolicy)

  # Merge any experimental toggles passed in dots
  knownExpToggles <- c("blockBanded", "factorSchurMultiFactor", "dFirstML",
                       "matrixFreePcg", "cholmodFinalSelected", "pdCertificates", "panelCi2", "profile", "blockSchur", "traceBatch", "sparseFactors", "saddlePoint")
  for(toggle in knownExpToggles){
    if(toggle %in% names(dots) && !toggle %in% names(experimental)){
      experimental[[toggle]] <- dots[[toggle]]
    }
  }

  # Experimental toggles
  expResolved <- list(
    blockBanded = .mme_experimental_mode(experimental$blockBanded, "auto", "sommer.mme.blockBanded"),
    factorSchurMultiFactor = .mme_experimental_mode(experimental$factorSchurMultiFactor, "auto", "sommer.mme.factorSchurMultiFactor"),
    dFirstML = .mme_experimental_mode(experimental$dFirstML, "auto", "sommer.mme.dFirstML"),
    matrixFreePcg = .mme_experimental_mode(experimental$matrixFreePcg, "auto", "sommer.mme.matrixFreePcg"),
    cholmodFinalSelected = .mme_experimental_mode(experimental$cholmodFinalSelected, "auto", "sommer.mme.cholmodFinalSelected"),
    pdCertificates = .mme_experimental_mode(experimental$pdCertificates, "auto", "sommer.mme.pdCertificates"),
    profile = .mme_experimental_mode(experimental$profile, "off", "sommer.mme.profile"),
    blockSchur = .mme_experimental_mode(experimental$blockSchur, "auto", "sommer.mme.blockSchur"),
    traceBatch = .mme_experimental_mode(experimental$traceBatch, "off", "sommer.mme.traceBatch"),
    sparseFactors = .mme_experimental_mode(experimental$sparseFactors, "off", "sommer.mme.sparseFactors"),
    saddlePoint = .mme_experimental_mode(experimental$saddlePoint, "off", "sommer.mme.saddlePoint"),
    panelCi2 = if(!is.null(experimental$panelCi2)) as.integer(experimental$panelCi2) else 64L
  )
  if(!is.finite(expResolved$panelCi2) || expResolved$panelCi2 < 1L){
    stop("panelCi2 must be a positive integer.", call.=FALSE)
  }

  # PCG parameters
  pcgDefaults <- list(
    tol = 1.0e-8,
    maxIters = 0L,
    traceProbes = 8L,
    lanczosSteps = 20L,
    preconditioner = "diagonal",
    nystromRank = 32L
  )
  for(pcgArg in names(pcgDefaults)){
    if(pcgArg %in% names(dots) && !pcgArg %in% names(pcg)){
      pcg[[pcgArg]] <- dots[[pcgArg]]
    }
  }
  pcgResolved <- .mme_pcg_control(pcg)

  control <- list(
    solver = solver,
    gpu = gpu,
    gpuMinSize = as.integer(gpuMinSize),
    denseMemoryMB = denseMemoryMB,
    memoryPolicy = memoryPolicy,
    experimental = expResolved,
    pcg = pcgResolved,
    acceleration = acceleration,
    factorScoreAugmentation = factorScoreAugmentation,
    recordDiagnostics = isTRUE(recordDiagnostics),
    verboseDiagnostics = isTRUE(verboseDiagnostics)
  )
  class(control) <- "sommer_mmes_control"
  control
}

#' Print method for sommer_mmes_control
#' @param x An object of class \code{sommer_mmes_control}.
#' @param ... Ignored.
#' @method print sommer_mmes_control
#' @export
print.sommer_mmes_control <- function(x, ...){
  cat("-- sommer MME Engine Control --\n")
  cat("Solver backend       :", x$solver, "\n")
  cat("GPU dense kernels    :", isTRUE(x$gpu), "\n")
  cat("Dense memory budget  :", x$denseMemoryMB, "MB\n")
  cat("Memory policy        :", x$memoryPolicy, "\n")
  cat("Acceleration         :", x$acceleration, "\n")
  cat("Factor augmentation  :", x$factorScoreAugmentation, "\n")
  cat("Iterative (PCG)      : tol =", x$pcg$tol,
      "| probes =", x$pcg$traceProbes,
      "| Lanczos =", x$pcg$lanczosSteps,
      "| precond =", x$pcg$preconditioner, "\n")
  cat("Experimental paths   :\n")
  for(toggle in names(x$experimental)){
    cat("  -", sprintf("%-22s", paste0(toggle, ":")), x$experimental[[toggle]], "\n")
  }
  invisible(x)
}

# Internal helper to merge engine argument with mmes top-level arguments
.merge_mmes_control <- function(engine = NULL,
                                solver = "auto",
                                ...){
  dots <- list(...)
  if(inherits(engine, "sommer_mmes_control")){
    ctrl <- engine
    previousPcg <- ctrl$pcg
    gpuControl <- mmes_control(
      gpu = if("gpu" %in% names(dots)) dots$gpu else ctrl$gpu,
      gpuMinSize = if("gpuMinSize" %in% names(dots)) dots$gpuMinSize else ctrl$gpuMinSize)
    ctrl$gpu <- gpuControl$gpu
    ctrl$gpuMinSize <- gpuControl$gpuMinSize
    if(!missing(solver) && !is.null(solver) && !identical(solver, "auto")){
      ctrl$solver <- tolower(solver[1L])
    }
    if("acceleration" %in% names(dots) && !is.null(dots$acceleration)){
      ctrl$acceleration <- match.arg(tolower(dots$acceleration[1L]), c("none", "aitken", "hybrid-ai"))
    }
    if("factorScoreAugmentation" %in% names(dots) && !is.null(dots$factorScoreAugmentation)){
      ctrl$factorScoreAugmentation <- match.arg(tolower(dots$factorScoreAugmentation[1L]),
                                                c("none", "fixed-shape", "profile"))
    }
    if("pcgTol" %in% names(dots) && !is.null(dots$pcgTol)){
      ctrl$pcg$tol <- dots$pcgTol
    }
    if("pcgMaxIters" %in% names(dots) && !is.null(dots$pcgMaxIters)){
      ctrl$pcg$maxIters <- dots$pcgMaxIters
    }
    if("pcgTraceProbes" %in% names(dots) && !is.null(dots$pcgTraceProbes)){
      ctrl$pcg$traceProbes <- dots$pcgTraceProbes
    }
    if("pcgLanczosSteps" %in% names(dots) && !is.null(dots$pcgLanczosSteps)){
      ctrl$pcg$lanczosSteps <- dots$pcgLanczosSteps
    }
    if("pcgPreconditioner" %in% names(dots) && !is.null(dots$pcgPreconditioner)){
      ctrl$pcg$preconditioner <- dots$pcgPreconditioner
    }
    if("pcgNystromRank" %in% names(dots) && !is.null(dots$pcgNystromRank)){
      ctrl$pcg$nystromRank <- dots$pcgNystromRank
    }
    if("pcg" %in% names(dots)){
      if(!is.list(dots$pcg) || (length(dots$pcg) &&
         (is.null(names(dots$pcg)) || anyNA(names(dots$pcg)) ||
          any(!nzchar(names(dots$pcg))) || anyDuplicated(names(dots$pcg))))){
        stop("pcg must be a list with unique, non-empty names.", call.=FALSE)
      }
      ctrl$pcg <- utils::modifyList(ctrl$pcg, dots$pcg, keep.null=TRUE)
    }
    for(pair in list(c("traceProbes", "maxTraceProbes"),
                     c("lanczosSteps", "maxLanczosSteps"))){
      if(!pair[2L] %in% names(dots$pcg) &&
         identical(previousPcg[[pair[1L]]], previousPcg[[pair[2L]]]) &&
         !identical(ctrl$pcg[[pair[1L]]], previousPcg[[pair[1L]]])){
        ctrl$pcg[pair[2L]] <- list(NULL)
      }
    }
    ctrl$pcg <- .mme_pcg_control(ctrl$pcg)
    if("denseMemoryMB" %in% names(dots) && !is.null(dots$denseMemoryMB)){
      ctrl$denseMemoryMB <- as.double(dots$denseMemoryMB)
    }
    if("memoryPolicy" %in% names(dots) && !is.null(dots$memoryPolicy)){
      ctrl$memoryPolicy <- mmes_control(memoryPolicy = dots$memoryPolicy)$memoryPolicy
    }
    if("experimental" %in% names(dots) && is.list(dots$experimental)){
      ctrl$experimental <- utils::modifyList(ctrl$experimental, dots$experimental)
    }
    if("saddlePoint" %in% names(dots) &&
       !"saddlePoint" %in% names(dots$experimental)){
      ctrl$experimental$saddlePoint <- dots$saddlePoint
    }
    ctrl$experimental$saddlePoint <- .mme_experimental_mode(
      ctrl$experimental$saddlePoint, "off", "sommer.mme.saddlePoint")
    return(ctrl)
  }

  # Build from solver and dots
  pcgArgs <- list()
  if("pcgTol" %in% names(dots) && !is.null(dots$pcgTol)) pcgArgs$tol <- dots$pcgTol
  if("pcgMaxIters" %in% names(dots) && !is.null(dots$pcgMaxIters)) pcgArgs$maxIters <- dots$pcgMaxIters
  if("pcgTraceProbes" %in% names(dots) && !is.null(dots$pcgTraceProbes)) pcgArgs$traceProbes <- dots$pcgTraceProbes
  if("pcgLanczosSteps" %in% names(dots) && !is.null(dots$pcgLanczosSteps)) pcgArgs$lanczosSteps <- dots$pcgLanczosSteps
  if("pcgPreconditioner" %in% names(dots) && !is.null(dots$pcgPreconditioner)) pcgArgs$preconditioner <- dots$pcgPreconditioner
  if("pcgNystromRank" %in% names(dots) && !is.null(dots$pcgNystromRank)) pcgArgs$nystromRank <- dots$pcgNystromRank
  if("pcg" %in% names(dots)){
    .mme_pcg_control(dots$pcg, checkLimits=FALSE)
    pcgArgs <- utils::modifyList(pcgArgs, dots$pcg, keep.null=TRUE)
  }

  accel <- if("acceleration" %in% names(dots) && !is.null(dots$acceleration)) dots$acceleration else "none"
  fsAug <- if("factorScoreAugmentation" %in% names(dots) && !is.null(dots$factorScoreAugmentation)) dots$factorScoreAugmentation else "none"
  denseMB <- if("denseMemoryMB" %in% names(dots)) dots$denseMemoryMB else NULL
  expList <- if("experimental" %in% names(dots) && is.list(dots$experimental)) dots$experimental else list()

  knownExpToggles <- c("blockBanded", "factorSchurMultiFactor", "dFirstML",
                       "matrixFreePcg", "cholmodFinalSelected", "pdCertificates", "panelCi2", "profile", "blockSchur", "traceBatch", "sparseFactors", "saddlePoint")
  for(toggle in knownExpToggles){
    if(toggle %in% names(dots) && !toggle %in% names(expList)){
      expList[[toggle]] <- dots[[toggle]]
    }
  }

  mmes_control(
    solver = solver,
    gpu = if("gpu" %in% names(dots)) dots$gpu else NULL,
    gpuMinSize = if("gpuMinSize" %in% names(dots)) dots$gpuMinSize else NULL,
    denseMemoryMB = denseMB,
    memoryPolicy = if("memoryPolicy" %in% names(dots)) dots$memoryPolicy else NULL,
    experimental = expList,
    pcg = pcgArgs,
    acceleration = accel,
    factorScoreAugmentation = fsAug
  )
}

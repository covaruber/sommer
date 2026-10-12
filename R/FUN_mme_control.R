# Structured solver engine configuration for mmes()
# Encapsulates solver selection, memory budgets, iterative (PCG) parameters,
# and experimental path toggles without polluting global options.

#' Mixed Model Equation Engine Control
#'
#' Builds a structured control object defining solver algorithms, memory budgets,
#' iterative settings, and experimental numerical paths for \code{\link{mmes}}.
#'
#' @param solver Backend solver for coefficient matrix factorization. One of \code{"auto"},
#'   \code{"ldlt"}, \code{"cholmod"}, or \code{"pcg"}.
#' @param denseMemoryMB Upper bound on dense memory allocations (in megabytes) used by
#'   Schur-complement engines, block-chain / banded paths, and selected-inverse caches. Default is 1000.
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
#'   }
#' @param pcg A named list of iterative PCG solver settings:
#'   \itemize{
#'     \item \code{tol}: Convergence tolerance for iterative CG solves (default 1e-8).
#'     \item \code{maxIters}: Maximum iterations for PCG solves (0 = automatic, default 0).
#'     \item \code{traceProbes}: Number of stochastic probes for Hutchinson traces (default 8).
#'     \item \code{lanczosSteps}: Lanczos quadrature steps for stochastic log determinant (default 20).
#'     \item \code{preconditioner}: Preconditioner type, \code{"diagonal"} or \code{"nystrom"}.
#'     \item \code{nystromRank}: Target rank for Nyström preconditioner (default 32).
#'   }
#' @param acceleration Convergence accelerator. One of \code{"none"} or \code{"aitken"}.
#' @param factorScoreAugmentation Factor-score augmentation strategy. One of \code{"none"},
#'   \code{"fixed-shape"}, or \code{"profile"}.
#' @param recordDiagnostics Logical; whether to capture detailed engine diagnostics in the fitted object. Default is TRUE.
#' @param verboseDiagnostics Logical; whether to print diagnostic details during execution. Default is FALSE.
#' @param ... Additional experimental toggles or legacy options passed as named arguments.
#'
#' @return An object of class \code{"sommer_mme_control"}.
#' @export
#' @aliases mme_control
mmes_control <- function(solver = c("auto", "ldlt", "cholmod", "pcg"),
                         denseMemoryMB = NULL,
                         experimental = list(),
                         pcg = list(),
                         acceleration = c("none", "aitken"),
                         factorScoreAugmentation = c("none", "fixed-shape", "profile"),
                         recordDiagnostics = TRUE,
                         verboseDiagnostics = FALSE,
                         ...){
  dots <- list(...)
  if(missing(solver) || is.null(solver)){
    solver <- "auto"
  } else {
    solver <- match.arg(tolower(solver[1L]), c("auto", "ldlt", "cholmod", "pcg"))
  }

  if(missing(acceleration) || is.null(acceleration)){
    acceleration <- "none"
  } else {
    acceleration <- match.arg(tolower(acceleration[1L]), c("none", "aitken"))
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

  # Helper to normalize three-state experimental modes ("auto", "off", "force")
  .normalize_mode <- function(val, default="auto", optName=NULL){
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

  # Merge any experimental toggles passed in dots
  knownExpToggles <- c("blockBanded", "factorSchurMultiFactor", "dFirstML",
                       "matrixFreePcg", "cholmodFinalSelected", "pdCertificates", "panelCi2")
  for(toggle in knownExpToggles){
    if(toggle %in% names(dots) && !toggle %in% names(experimental)){
      experimental[[toggle]] <- dots[[toggle]]
    }
  }

  # Experimental toggles
  expResolved <- list(
    blockBanded = .normalize_mode(experimental$blockBanded, "auto", "sommer.mme.blockBanded"),
    factorSchurMultiFactor = .normalize_mode(experimental$factorSchurMultiFactor, "auto", "sommer.mme.factorSchurMultiFactor"),
    dFirstML = .normalize_mode(experimental$dFirstML, "auto", "sommer.mme.dFirstML"),
    matrixFreePcg = .normalize_mode(experimental$matrixFreePcg, "auto", "sommer.mme.matrixFreePcg"),
    cholmodFinalSelected = .normalize_mode(experimental$cholmodFinalSelected, "auto", "sommer.mme.cholmodFinalSelected"),
    pdCertificates = .normalize_mode(experimental$pdCertificates, "auto", "sommer.mme.pdCertificates"),
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
  pcgResolved <- utils::modifyList(pcgDefaults, pcg)
  if(!is.numeric(pcgResolved$tol) || pcgResolved$tol <= 0){
    stop("pcg$tol must be a positive number.", call.=FALSE)
  }
  pcgResolved$preconditioner <- match.arg(tolower(pcgResolved$preconditioner[1L]),
                                          c("diagonal", "nystrom"))
  pcgResolved$maxIters <- as.integer(pcgResolved$maxIters)
  pcgResolved$traceProbes <- as.integer(pcgResolved$traceProbes)
  pcgResolved$lanczosSteps <- as.integer(pcgResolved$lanczosSteps)
  pcgResolved$nystromRank <- as.integer(pcgResolved$nystromRank)

  control <- list(
    solver = solver,
    denseMemoryMB = denseMemoryMB,
    experimental = expResolved,
    pcg = pcgResolved,
    acceleration = acceleration,
    factorScoreAugmentation = factorScoreAugmentation,
    recordDiagnostics = isTRUE(recordDiagnostics),
    verboseDiagnostics = isTRUE(verboseDiagnostics)
  )
  class(control) <- "sommer_mme_control"
  control
}

#' @rdname mmes_control
#' @export
mme_control <- mmes_control

#' Print method for sommer_mme_control
#' @param x An object of class \code{sommer_mme_control}.
#' @param ... Ignored.
#' @export
print.sommer_mme_control <- function(x, ...){
  cat("-- sommer MME Engine Control --\n")
  cat("Solver backend       :", x$solver, "\n")
  cat("Dense memory budget  :", x$denseMemoryMB, "MB\n")
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
  if(inherits(engine, "sommer_mme_control")){
    ctrl <- engine
    if(!missing(solver) && !is.null(solver) && !identical(solver, "auto")){
      ctrl$solver <- tolower(solver[1L])
    }
    if("acceleration" %in% names(dots) && !is.null(dots$acceleration)){
      ctrl$acceleration <- match.arg(tolower(dots$acceleration[1L]), c("none", "aitken"))
    }
    if("factorScoreAugmentation" %in% names(dots) && !is.null(dots$factorScoreAugmentation)){
      ctrl$factorScoreAugmentation <- match.arg(tolower(dots$factorScoreAugmentation[1L]),
                                                c("none", "fixed-shape", "profile"))
    }
    if("pcgTol" %in% names(dots) && !is.null(dots$pcgTol)){
      ctrl$pcg$tol <- as.double(dots$pcgTol)
    }
    if("pcgMaxIters" %in% names(dots) && !is.null(dots$pcgMaxIters)){
      ctrl$pcg$maxIters <- as.integer(dots$pcgMaxIters)
    }
    if("pcgTraceProbes" %in% names(dots) && !is.null(dots$pcgTraceProbes)){
      ctrl$pcg$traceProbes <- as.integer(dots$pcgTraceProbes)
    }
    if("pcgLanczosSteps" %in% names(dots) && !is.null(dots$pcgLanczosSteps)){
      ctrl$pcg$lanczosSteps <- as.integer(dots$pcgLanczosSteps)
    }
    if("pcgPreconditioner" %in% names(dots) && !is.null(dots$pcgPreconditioner)){
      ctrl$pcg$preconditioner <- match.arg(tolower(dots$pcgPreconditioner[1L]), c("diagonal", "nystrom"))
    }
    if("pcgNystromRank" %in% names(dots) && !is.null(dots$pcgNystromRank)){
      ctrl$pcg$nystromRank <- as.integer(dots$pcgNystromRank)
    }
    if("denseMemoryMB" %in% names(dots) && !is.null(dots$denseMemoryMB)){
      ctrl$denseMemoryMB <- as.double(dots$denseMemoryMB)
    }
    if("experimental" %in% names(dots) && is.list(dots$experimental)){
      ctrl$experimental <- utils::modifyList(ctrl$experimental, dots$experimental)
    }
    return(ctrl)
  }

  # Build from solver and dots
  pcgArgs <- list()
  if("pcgTol" %in% names(dots) && !is.null(dots$pcgTol)) pcgArgs$tol <- as.double(dots$pcgTol)
  if("pcgMaxIters" %in% names(dots) && !is.null(dots$pcgMaxIters)) pcgArgs$maxIters <- as.integer(dots$pcgMaxIters)
  if("pcgTraceProbes" %in% names(dots) && !is.null(dots$pcgTraceProbes)) pcgArgs$traceProbes <- as.integer(dots$pcgTraceProbes)
  if("pcgLanczosSteps" %in% names(dots) && !is.null(dots$pcgLanczosSteps)) pcgArgs$lanczosSteps <- as.integer(dots$pcgLanczosSteps)
  if("pcgPreconditioner" %in% names(dots) && !is.null(dots$pcgPreconditioner)) pcgArgs$preconditioner <- dots$pcgPreconditioner
  if("pcgNystromRank" %in% names(dots) && !is.null(dots$pcgNystromRank)) pcgArgs$nystromRank <- as.integer(dots$pcgNystromRank)
  if("pcg" %in% names(dots) && is.list(dots$pcg)) pcgArgs <- utils::modifyList(pcgArgs, dots$pcg)

  accel <- if("acceleration" %in% names(dots) && !is.null(dots$acceleration)) dots$acceleration else "none"
  fsAug <- if("factorScoreAugmentation" %in% names(dots) && !is.null(dots$factorScoreAugmentation)) dots$factorScoreAugmentation else "none"
  denseMB <- if("denseMemoryMB" %in% names(dots)) dots$denseMemoryMB else NULL
  expList <- if("experimental" %in% names(dots) && is.list(dots$experimental)) dots$experimental else list()

  knownExpToggles <- c("blockBanded", "factorSchurMultiFactor", "dFirstML",
                       "matrixFreePcg", "cholmodFinalSelected", "pdCertificates", "panelCi2")
  for(toggle in knownExpToggles){
    if(toggle %in% names(dots) && !toggle %in% names(expList)){
      expList[[toggle]] <- dots[[toggle]]
    }
  }

  mmes_control(
    solver = solver,
    denseMemoryMB = denseMB,
    experimental = expList,
    pcg = pcgArgs,
    acceleration = accel,
    factorScoreAugmentation = fsAug
  )
}

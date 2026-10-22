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

.mme_compression_control <- function(value){
  if(!is.character(value) || length(value) != 1L || is.na(value) ||
     !tolower(value) %in% c("off", "groups", "classes")){
    stop("representation$compression must be one of 'off', 'groups', or 'classes'.", call.=FALSE)
  }
  tolower(value)
}

.mme_covariance_information <- function(value){
  if(!is.character(value) || length(value) != 1L || is.na(value) ||
     !tolower(value) %in% c("full", "none")){
    stop("output$covariance must be 'full' or 'none'.", call.=FALSE)
  }
  tolower(value)
}

.mme_pcg_mapping <- function(){
  c(tol="solve$tol", maxIters="solve$maxIters",
    preconditioner="solve$preconditioner", nystromRank="solve$nystromRank",
    warmStart="solve$warmStart", recycleRank="solve$recycleRank",
    recycleMethod="solve$recycleMethod", recycleState="solve$recycleState",
    exportRecycle="solve$exportRecycle", rhsMethod="solve$rhsMethod",
    blockSize="solve$blockSize", residualRefresh="solve$residualRefresh",
    logDetScaling="logdet$scaling", lanczosSteps="logdet$lanczosSteps",
    traceProbes="score$traceProbes", traceEstimator="score$estimator",
    traceRank="score$rank", accuracy="convergence$accuracy",
    maxTraceProbes="convergence$maxTraceProbes",
    maxLanczosSteps="convergence$maxLanczosSteps", logDetTol="convergence$logDetTol",
    scoreTol="convergence$scoreTol", seed="execution$seed",
    compression="representation$compression")
}

.mme_pcg_public <- function(name){
  paths <- .mme_pcg_mapping()
  if(!name %in% names(paths)){
    stop("Unknown internal PCG setting: ", name, ".", call.=FALSE)
  }
  paths[[name]]
}

.mme_pcg_control <- function(pcg = list(), checkLimits = TRUE){
  defaults <- list(
    tol = 1e-8, maxIters = 0L, traceProbes = 8L, lanczosSteps = 20L,
    preconditioner = "diagonal", nystromRank = 32L, warmStart = FALSE,
    recycleRank = 0L, recycleMethod = "coarse",
    recycleState = NULL, exportRecycle = FALSE, compression = "off",
    rhsMethod = "independent", blockSize = 16L,
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
      stop(.mme_pcg_public(name), " must be one positive finite number.", call.=FALSE)
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
      stop(.mme_pcg_public(name), " must be one integer between ", minima[[name]],
           " and ", .Machine$integer.max, ".", call.=FALSE)
    }
    out[[name]] <- as.integer(value)
  }
  choices <- list(
    preconditioner=c("diagonal", "nystrom", "block-jacobi"),
    recycleMethod=c("coarse", "krylov"),
    compression=c("off", "groups", "classes"),
    rhsMethod=c("independent", "block"),
    logDetScaling=c("none", "diagonal"),
    traceEstimator=c("hutchinson", "control-variate"),
    accuracy=c("fixed", "adaptive")
  )
  for(name in names(choices)){
    value <- out[[name]]
    if(!is.character(value) || length(value) != 1L || is.na(value) ||
       !tolower(value) %in% choices[[name]]){
      stop(.mme_pcg_public(name), " must be one of ",
           paste(shQuote(choices[[name]]), collapse=", "), ".", call.=FALSE)
    }
    out[[name]] <- tolower(value)
  }
  if(!is.logical(out$warmStart) || length(out$warmStart) != 1L ||
     is.na(out$warmStart)){
    stop("solve$warmStart must be a single TRUE/FALSE value.", call.=FALSE)
  }
  if(!is.logical(out$exportRecycle) || length(out$exportRecycle) != 1L ||
     is.na(out$exportRecycle)){
    stop("solve$exportRecycle must be a single TRUE/FALSE value.", call.=FALSE)
  }
  if(!is.null(out$recycleState)){
    state <- out$recycleState
    if(!inherits(state, "sommer_pcg_recycle") || !is.list(state) ||
       !identical(names(state), c("version", "coordinates", "basis")) ||
       !identical(state$version, 1L) || !is.list(state$coordinates) ||
       !is.matrix(state$basis) || !is.numeric(state$basis) ||
       is.complex(state$basis) || !nrow(state$basis) || !ncol(state$basis) ||
       any(!is.finite(state$basis))){
      stop("solve$recycleState must be a valid exported sommer_pcg_recycle state.",
           call.=FALSE)
    }
    if(ncol(state$basis) > out$recycleRank){
      stop("solve$recycleState rank exceeds solve$recycleRank.", call.=FALSE)
    }
  }
  if((out$exportRecycle || !is.null(out$recycleState)) &&
     (out$recycleMethod != "krylov" || out$recycleRank == 0L)){
    stop("solve$recycleState and solve$exportRecycle require solve$recycleMethod='krylov' and positive solve$recycleRank.",
         call.=FALSE)
  }
  if(checkLimits && (out$maxTraceProbes < out$traceProbes ||
     out$maxLanczosSteps < out$lanczosSteps)){
    stop("convergence$maxTraceProbes and convergence$maxLanczosSteps must not be smaller than score$traceProbes and logdet$lanczosSteps.",
         call.=FALSE)
  }
  out
}

# Private validated transport for the native descriptor contract.
# Public method ownership and strict module validation live in FUN_mme_modules.R.
.mmes_engine_control <- function(solver = c("auto", "ldlt", "cholmod", "pcg"),
                         denseMemoryMB = NULL,
                         experimental = list(),
                         pcg = list(),
                         acceleration = c("none", "aitken", "hybrid-ai"),
                         factorScoreAugmentation = c("none", "fixed-shape", "profile"),
                         memoryPolicy = NULL,
                         gpu = NULL,
                         gpuMinSize = NULL,
                         compression = "off",
                         covarianceInformation = "full",
                         retainCoefficientMatrix = TRUE){
  if(!is.logical(retainCoefficientMatrix) || length(retainCoefficientMatrix) != 1L ||
     is.na(retainCoefficientMatrix)){
    stop("retainCoefficientMatrix must be a single TRUE/FALSE value.", call.=FALSE)
  }
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
    sparseFactors = .mme_experimental_mode(experimental$sparseFactors, "auto", "sommer.mme.sparseFactors"),
    saddlePoint = .mme_experimental_mode(experimental$saddlePoint, "off", "sommer.mme.saddlePoint"),
    reducedAI = .mme_experimental_mode(experimental$reducedAI, "off", "sommer.mme.reducedAI"),
    commonScale = .mme_experimental_mode(experimental$commonScale, "off", "sommer.mme.commonScale"),
    panelCi2 = if(!is.null(experimental$panelCi2)) as.integer(experimental$panelCi2) else 64L
  )
  if(!is.finite(expResolved$panelCi2) || expResolved$panelCi2 < 1L){
    stop("panelCi2 must be a positive integer.", call.=FALSE)
  }

  pcgResolved <- .mme_pcg_control(pcg)
  compression <- .mme_compression_control(compression)
  covarianceInformation <- .mme_covariance_information(covarianceInformation)

  control <- list(
    solver = solver,
    gpu = gpu,
    gpuMinSize = as.integer(gpuMinSize),
    denseMemoryMB = denseMemoryMB,
    memoryPolicy = memoryPolicy,
    experimental = expResolved,
    pcg = pcgResolved,
    compression = compression,
    covarianceInformation = covarianceInformation,
    retainCoefficientMatrix = retainCoefficientMatrix,
    acceleration = acceleration,
    factorScoreAugmentation = factorScoreAugmentation
  )
  class(control) <- "sommer_mmes_engine"
  control
}

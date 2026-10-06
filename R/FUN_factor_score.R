.factor_score_decomposition <- function(factor){
  q <- as.integer(factor$dim)
  par <- as.numeric(factor$par)
  representation <- factor$representation$kind

  if(identical(representation, "lowrank_diagonal") &&
     factor$model %in% c("fa", "rr")){
    order <- as.integer(factor$order)
    nload <- if(factor$model == "fa") as.integer(factor$fa_nload) else
      as.integer(factor$rr_nload)
    rows <- if(factor$model == "fa") as.integer(factor$fa_row) else
      as.integer(factor$rr_row)
    cols <- if(factor$model == "fa") as.integer(factor$fa_col) else
      as.integer(factor$rr_col)
    diagonal <- if(factor$model == "fa") as.logical(factor$fa_diag) else
      as.logical(factor$rr_diag)
    loading <- matrix(0, q, order)
    for(index in seq_len(nload)){
      loading[rows[index], cols[index]] <-
        if(diagonal[index]) exp(par[index]) else par[index]
    }
    specific <- if(factor$model == "fa"){
      c(1, exp(par[nload + seq_len(q - 1L)]))
    }else{
      rep(1, q)
    }
    return(list(loading=loading, specific=specific,
                normalization=1 + loading[1L, 1L]^2))
  }

  if(identical(representation, "compound_symmetry")){
    lower <- factor$report$lower[1L]
    upper <- factor$report$upper[1L]
    rho <- lower + (upper - lower) * stats::plogis(par[1L])
    if(!length(par) %in% c(1L, q) || !is.finite(rho) || rho < 0 || rho >= 1){
      return(NULL)
    }
    variance <- if(length(par) == 1L) rep(1, q) else c(1, exp(par[-1L]))
    loading <- matrix(sqrt(rho * variance), ncol=1L)
    return(list(loading=loading, specific=variance * (1 - rho),
                normalization=1))
  }

  NULL
}

.factor_score_augment <- function(Z, Ai, covStruct, Zind, rtermss, rTermsNames,
                                 allowFree=FALSE){
  expandedZ <- list()
  expandedAi <- list()
  expandedCovStruct <- list()
  expandedZind <- numeric()
  expandedTerms <- character()
  expandedTermNames <- list()
  mappings <- list()
  constrainedStructures <- list()
  oldTerms <- length(covStruct)

  for(term in seq_len(oldTerms)){
    zIndex <- which(Zind == term)
    shape <- covStruct[[term]]
    factors <- shape$factors
    representation <- if(length(factors) == 1L){
      factors[[1L]]$representation$kind
    }else{
      NULL
    }
    if(is.null(representation) ||
       !representation %in% c("lowrank_diagonal", "compound_symmetry")){
      expandedZ <- c(expandedZ, Z[zIndex])
      newIndex <- length(expandedCovStruct) + 1L
      expandedAi[[newIndex]] <- Ai[[term]]
      expandedCovStruct[[newIndex]] <- shape
      expandedZind <- c(expandedZind, rep(newIndex, length(zIndex)))
      expandedTerms <- c(expandedTerms, rtermss[term])
      expandedTermNames[[newIndex]] <- rTermsNames[[term]]
      mappings[[term]] <- list(newIndices=newIndex, augmented=FALSE,
               q=as.integer(shape$dim), d=ncol(Z[[zIndex[1L]]]),
               levels=as.character(shape$levels))
      next
    }

    factor <- factors[[1L]]
    factorStart <- as.integer(factor$par_start)
    factorEnd <- as.integer(factor$par_end)
    factor$par <- if(factorEnd >= factorStart) shape$par[factorStart:factorEnd] else numeric()
    if(any(shape$free[-1L]) && !allowFree){
      stop("factorScoreAugmentation='fixed-shape' requires all covariance-shape parameters to be fixed; only the shared variance scale may be estimated.", call.=FALSE)
    }
    q <- as.integer(factor$dim)
    if(length(zIndex) != q || any(vapply(Z[zIndex], ncol, integer(1)) != ncol(Z[[zIndex[1L]]]))){
      stop("Factor-score augmentation requires one equal-width incidence block per covariance level.", call.=FALSE)
    }
    decomposition <- .factor_score_decomposition(factor)
    if(is.null(decomposition)){
      stop("Factor-score augmentation requires a positive low-rank-plus-diagonal covariance representation.", call.=FALSE)
    }
    loading <- decomposition$loading
    specific <- decomposition$specific
    normalization <- decomposition$normalization
    order <- ncol(loading)
    sourceDesigns <- Z[zIndex]
    latentDesigns <- lapply(seq_len(order), function(component){
      Reduce(`+`, Map(function(design, coefficient) design * coefficient,
                      sourceDesigns, loading[, component]))
    })
    specificDesigns <- Map(function(design, variance) design * sqrt(variance),
                           sourceDesigns, specific)
    designs <- c(latentDesigns, specificDesigns)
    shapeScale <- shape$par[1L] - log(normalization)
    newIndices <- integer(length(designs))
    for(component in seq_along(designs)){
      newIndex <- length(expandedCovStruct) + 1L
      newIndices[component] <- newIndex
      expandedZ[[length(expandedZ) + 1L]] <- designs[[component]]
      expandedAi[[newIndex]] <- Ai[[term]]
      label <- if(component <= order) paste0("factor", component) else
        paste0("specific", component - order)
      expandedCovStruct[[newIndex]] <- list(
        type="kron", descriptor_version=2L, dim=1L, levels=label,
        par=stats::setNames(shapeScale, "sigma2"), free=TRUE,
        par_names=paste0(rtermss[term], ":", label, ":sigma2"),
        factors=list(), sigma2_is_default=FALSE
      )
      expandedZind <- c(expandedZind, newIndex)
      expandedTerms <- c(expandedTerms, paste0(rtermss[term], "[", label, "]"))
      expandedTermNames[[newIndex]] <- paste0(rtermss[term], ":", label, ":sigma2")
    }
    constrainedStructures[[length(constrainedStructures) + 1L]] <- newIndices
    mappings[[term]] <- list(newIndices=newIndices, augmented=TRUE, loading=loading,
                             specific=specific, normalization=normalization,
                             levels=as.character(shape$levels), d=ncol(sourceDesigns[[1L]]),
                             shape=shape, sourceDesigns=sourceDesigns,
                             originalTerm=rtermss[term])
  }

  list(Z=expandedZ, Ai=expandedAi, covStruct=expandedCovStruct,
       Zind=expandedZind, terms=expandedTerms, termNames=expandedTermNames,
       mappings=mappings, constrainedStructures=constrainedStructures,
      originalTerms=rtermss, originalTermNames=rTermsNames,
      originalCovStruct=covStruct)
}

.factor_score_restore <- function(result, info, X, originalZ, residualCovStruct){
  result$C_augmented <- result$C
  result$bu_augmented <- result$bu
  result$uListAugmented <- result$uList
  result$partitionsAugmented <- result$partitions

  originalU <- vector("list", length(info$mappings))
  originalTheta <- vector("list", length(info$mappings) + 1L)
  originalCovPar <- vector("list", length(info$mappings) + 1L)
  originalCovStruct <- info$originalCovStruct
  originalPartitions <- vector("list", length(info$mappings))
  monitorRows <- vector("list", length(info$mappings))
  vcStart <- cumsum(c(1L, vapply(result$covPar, length, integer(1))))
  publicStart <- ncol(X) + 1L
  augmentedCoefficientMap <- vector("list", length(info$mappings))

  for(term in seq_along(info$mappings)){
    mapping <- info$mappings[[term]]
    newIndices <- mapping$newIndices
    shape <- info$originalCovStruct[[term]]
    q <- as.integer(shape$dim)
    d <- mapping$d
    partition <- cbind(
      publicStart + (seq_len(q) - 1L) * d,
      publicStart + seq_len(q) * d - 1L
    )
    originalPartitions[[term]] <- partition
    mapping$publicStart <- publicStart
    mapping$publicRanges <- partition
    mapping$augmentedRanges <- lapply(newIndices, function(index){
      result$partitionsAugmented[[index]]
    })
    augmentedCoefficientMap[[term]] <- mapping

    if(mapping$augmented){
      factorCount <- ncol(mapping$loading)
      factorScores <- matrix(unlist(lapply(newIndices[seq_len(factorCount)],
        function(index) result$uListAugmented[[index]][, 1L]), use.names=FALSE),
        nrow=d, ncol=factorCount)
      specificScores <- matrix(unlist(lapply(newIndices[factorCount + seq_len(q)],
        function(index) result$uListAugmented[[index]][, 1L]), use.names=FALSE),
        nrow=d, ncol=q)
      originalU[[term]] <- factorScores %*% t(mapping$loading) +
        specificScores %*% diag(sqrt(mapping$specific), nrow=q)

      sigmaAugmented <- as.numeric(result$theta[[newIndices[1L]]][1L,1L])
      loading <- mapping$loading
      covariance <- tcrossprod(loading) + diag(mapping$specific, q)
      originalTheta[[term]] <- sigmaAugmented * covariance
      sigmaWorking <- result$covParWorking[[newIndices[1L]]][1L] +
        log(mapping$normalization)
      originalCovStruct[[term]]$par[1L] <- sigmaWorking
      originalCovStruct[[term]]$sigma2_is_default <- FALSE
      sigmaOriginal <- sigmaAugmented * mapping$normalization
      factorPar <- shape$par[-1L]
      for(position in seq_along(factorPar)){
        kind <- .vc_kind(shape, position + 1L)
        factorPar[position] <- .vc_to_natural(factorPar[position], kind$kind,
                                               kind$lower, kind$upper)
      }
      originalCovPar[[term]] <- c(sigmaOriginal,
        stats::setNames(factorPar, shape$par_names[-1L]))
      monitorBase <- result$monitor[vcStart[newIndices[1L]], , drop=FALSE]
      monitorRows[[term]] <- rbind(
        as.numeric(monitorBase) + log(mapping$normalization),
        matrix(shape$par[-1L], nrow=length(shape$par)-1L,
               ncol=ncol(monitorBase))
      )
      originalCovStruct[[term]]$par[-1L] <- shape$par[-1L]
      originalCovStruct[[term]]$par[1L] <- log(sigmaOriginal)
    }else{
      originalU[[term]] <- result$uListAugmented[[newIndices[1L]]]
      originalTheta[[term]] <- result$theta[[newIndices[1L]]]
      originalCovPar[[term]] <- result$covPar[[newIndices[1L]]]
      monitorRows[[term]] <- result$monitor[
        vcStart[newIndices[1L]]:(vcStart[newIndices[1L] + 1L] - 1L), , drop=FALSE]
    }
    publicStart <- publicStart + q * d
  }

  residualIndex <- length(result$theta)
  residualFirstVC <- vcStart[residualIndex]
  originalMonitor <- do.call(rbind, c(monitorRows,
    list(result$monitor[residualFirstVC:(nrow(result$monitor)), , drop=FALSE])))
  originalTheta[[length(originalTheta)]] <- result$theta[[residualIndex]]
  originalCovPar[[length(originalCovPar)]] <- result$covPar[[residualIndex]]
  originalCovStruct[[length(originalCovStruct) + 1L]] <- residualCovStruct

  originalSECount <- sum(vapply(originalCovPar, length, integer(1)))
  originalSEMap <- matrix(0, originalSECount, nrow(result$theta_se))
  publicOffset <- 0L
  for(term in seq_along(info$mappings)){
    mapping <- info$mappings[[term]]
    newIndex <- mapping$newIndices
    if(mapping$augmented){
      sigmaAugmented <- as.numeric(result$theta[[newIndex[1L]]][1L,1L])
      originalSEMap[publicOffset + 1L, vcStart[newIndex[1L]]] <-
        mapping$normalization
    }else{
      npar <- length(originalCovPar[[term]])
      oldStart <- vcStart[newIndex[1L]]
      originalSEMap[publicOffset + seq_len(npar), oldStart:(oldStart+npar-1L)] <- diag(npar)
    }
    publicOffset <- publicOffset + length(originalCovPar[[term]])
  }
  residualPar <- length(originalCovPar[[length(originalCovPar)]])
  originalSEMap[publicOffset + seq_len(residualPar),
                residualFirstVC:(residualFirstVC+residualPar-1L)] <- diag(residualPar)

  originalUVector <- unlist(lapply(originalU, as.vector), use.names=FALSE)
  originalNames <- unlist(lapply(originalZ, colnames), use.names=FALSE)
  result$uList <- originalU
  result$u <- matrix(originalUVector, ncol=1L,
                     dimnames=list(originalNames, NULL))
  result$bu <- rbind(result$b, result$u)
  result$W <- do.call(cbind, c(list(X), originalZ))
  result$partitions <- originalPartitions
  result$theta <- originalTheta
  result$covPar <- originalCovPar
  result$covStruct <- originalCovStruct
  result$monitor <- originalMonitor
  result$theta_se <- originalSEMap %*% result$theta_se %*% t(originalSEMap)
  result$factorScoreInfo <- info
  result$factorScoreInfo$mappings <- augmentedCoefficientMap
  result$CRepresentation <- "factor-score-augmented"
  result
}

.mmes_factor_score_profile <- function(call, caller, family, nIters, REML,
                                       henderson, computeCi, vcc, solver){
  if(!inherits(family, "family") || !identical(family$family, "gaussian") ||
     !identical(family$link, "identity") || !isTRUE(REML) || !isTRUE(henderson) ||
      computeCi != 0L || !is.null(vcc) ||
      !(tolower(solver) %in% c("auto", "ldlt", "cholmod"))){
    stop("factorScoreAugmentation='profile' currently requires Gaussian REML Henderson fitting, computeCi=0, and no user vcc constraints or rotation.", call.=FALSE)
  }
  setupCall <- call
  setupCall$factorScoreAugmentation <- "none"
  setupCall$returnParam <- TRUE
  setup <- eval(setupCall, envir=caller)
  if(!is.null(setup$rotation)){
    stop("factorScoreAugmentation='profile' does not support rotation=TRUE.", call.=FALSE)
  }
  shapeTerms <- integer()
  startingShapes <- list()
  freeShapes <- list()
  for(term in seq_along(setup$rtermss)){
    shape <- setup$covStruct[[term]]
    if(length(shape$factors) != 1L ||
       !shape$factors[[1L]]$representation$kind %in%
         c("lowrank_diagonal", "compound_symmetry")) next
    free <- which(shape$free[-1L])
    if(!length(free)) next
    shapeTerms <- c(shapeTerms, term)
    startingShapes[[as.character(term)]] <- shape$par[-1L]
    freeShapes[[as.character(term)]] <- free
  }
  if(!length(shapeTerms)){
    stop("factorScoreAugmentation='profile' requires free shape parameters on an eligible low-rank or compound-symmetry random term.", call.=FALSE)
  }

  initial <- unlist(lapply(names(startingShapes), function(term){
    startingShapes[[term]][freeShapes[[term]]]
  }), use.names=FALSE)

  unpack <- function(parameters){
    updated <- startingShapes
    offset <- 0L
    for(term in names(updated)){
      positions <- freeShapes[[term]]
      count <- length(positions)
      updated[[term]][positions] <- parameters[offset + seq_len(count)]
      offset <- offset + count
    }
    out <- vector("list", max(shapeTerms))
    for(term in shapeTerms) out[[term]] <- updated[[as.character(term)]]
    out
  }
  augmentedSetupCall <- call
  augmentedSetupCall$factorScoreAugmentation <- "fixed-shape"
  augmentedSetupCall$.factorScoreParameters <- as.name(".factorScoreParameters")
  augmentedSetupCall$returnParam <- TRUE
  augmentedSetupEnvironment <- new.env(parent=caller)
  augmentedSetupEnvironment$.factorScoreParameters <- unpack(initial)
  augmentedSetup <- eval(augmentedSetupCall, envir=augmentedSetupEnvironment)
  for(group in augmentedSetup$factorScoreInfo$constrainedStructures){
    for(structure in group){
      augmentedSetup$covStruct[[structure]]$free[1L] <- TRUE
    }
  }
  innerStart <- list(parameters=lapply(augmentedSetup$covStruct, `[[`, "par"),
    covStruct=augmentedSetup$covStruct, included=augmentedSetup$obsInfo$included)
  cachedParameters <- NULL
  cachedObjective <- NULL
  bestParameters <- NULL
  bestValue <- Inf
  bestFit <- NULL
  profileEvaluations <- 0L
  validEvaluations <- 0L
  invalidEvaluations <- 0L
  totalInnerIterations <- 0L
  totalSymbolicAnalyses <- 0L
  warmStartedEvaluations <- 0L
  objective <- function(parameters){
    if(!is.null(cachedParameters) && identical(as.numeric(parameters), cachedParameters)){
      return(cachedObjective)
    }
    innerCall <- call
    innerCall$factorScoreAugmentation <- "fixed-shape"
    innerCall$.factorScoreParameters <- as.name(".factorScoreParameters")
    innerCall$.pqlInner <- TRUE
    innerCall$.pqlStart <- as.name(".factorScoreStart")
    profileEnvironment <- new.env(parent=caller)
    profileEnvironment$.factorScoreParameters <- unpack(parameters)
    profileEnvironment$.factorScoreStart <- innerStart
    profileEvaluations <<- profileEvaluations + 1L
    fit <- tryCatch(eval(innerCall, envir=profileEnvironment), error=function(e) NULL)
    if(is.null(fit) || !length(fit$llik) || !is.finite(tail(as.numeric(fit$llik),1L))){
      invalidEvaluations <<- invalidEvaluations + 1L
      return(1e100)
    }
    validEvaluations <<- validEvaluations + 1L
    totalInnerIterations <<- totalInnerIterations + ncol(fit$monitor)
    if(isTRUE(fit$pqlWarmStarted) || isTRUE(fit$factorProfileWarmStarted)){
      warmStartedEvaluations <<- warmStartedEvaluations + 1L
    }
    if(!is.null(fit$engineDiagnostics$CsymbolicAnalyses)){
      totalSymbolicAnalyses <<- totalSymbolicAnalyses + fit$engineDiagnostics$CsymbolicAnalyses
    }
    cachedParameters <<- as.numeric(parameters)
    cachedObjective <<- -tail(as.numeric(fit$llik),1L)
    if(cachedObjective < bestValue){
      bestValue <<- cachedObjective
      bestParameters <<- cachedParameters
      bestFit <<- fit
    }
    if(!is.null(fit$covParWorking)){
      innerStart <<- list(parameters=fit$covParWorking,
        covStruct=augmentedSetup$covStruct, included=fit$obsInfo$included,
        ldltCache=fit$.ldltCache, cholmodCache=fit$.cholmodCache)
    }
    cachedObjective
  }
  lower <- upper <- rep(6, length(initial))
  for(term in names(freeShapes)){
    shape <- setup$covStruct[[as.integer(term)]]
    factor <- shape$factors[[1L]]
    positions <- freeShapes[[term]]
    offset <- sum(vapply(names(freeShapes)[seq_len(match(term,names(freeShapes))-1L)],
      function(previous) length(freeShapes[[previous]]), integer(1)))
    if(identical(factor$representation$kind, "compound_symmetry")){
      lower[offset + seq_along(positions)] <- -10
      upper[offset + seq_along(positions)] <- 10
      rhoPosition <- which(positions == 1L)
      if(length(rhoPosition)){
        rhoLower <- factor$report$lower[1L]
        rhoUpper <- factor$report$upper[1L]
        zeroRhoEta <- stats::qlogis((0 - rhoLower) / (rhoUpper - rhoLower))
        lower[offset + rhoPosition] <- max(lower[offset + rhoPosition], zeroRhoEta)
      }
    }else{
      nload <- if(factor$model == "fa") factor$fa_nload else factor$rr_nload
      for(local in seq_along(positions)){
        position <- positions[local]
        upper[offset + local] <- if(position <= nload) 6 else 10
        lower[offset + local] <- -upper[offset + local]
      }
    }
  }
  lower[1L] <- max(lower[1L], -6)
  startValue <- objective(initial)
  outerMaxIterations <- max(5L, min(12L, as.integer(600L / length(initial))))
  optimum <- stats::optim(initial, objective, method="L-BFGS-B",
    lower=lower, upper=upper,
    control=list(maxit=outerMaxIterations,
                 factr=1e7, pgtol=1e-5))
  if(is.null(bestFit)){
    stop("Factor-score profile optimizer failed to produce a valid augmented REML fit.", call.=FALSE)
  }
  result <- bestFit
  for(term in shapeTerms){
    original <- setup$covStruct[[term]]
    result$covStruct[[term]]$free <- original$free
    result$covStruct[[term]]$factors[[1L]]$free <- original$factors[[1L]]$free
  }
  result$vcParams <- .vc_param_table(result$covStruct,
    c(names(result$partitions), names(result$theta)[length(result$theta)]))
  result$factorScoreAugmentation <- "profile"
  result$factorScoreProfile <- list(convergence=optimum$convergence,
    message=optimum$message, objective=bestValue, startObjective=startValue,
    evaluations=optimum$counts, parameters=bestParameters,
    outerMaxIterations=outerMaxIterations,
    profileEvaluations=profileEvaluations, validEvaluations=validEvaluations,
    invalidEvaluations=invalidEvaluations, totalInnerIterations=totalInnerIterations,
    totalSymbolicAnalyses=totalSymbolicAnalyses,
    warmStartedEvaluations=warmStartedEvaluations,
    method="bounded L-BFGS-B outer profile; augmented REML inner fits")
  result$.ldltCache <- NULL
  result$.cholmodCache <- NULL
  attr(result$covStruct, "ldltCache") <- NULL
  attr(result$covStruct, "cholmodCache") <- NULL
  shapeRows <- unlist(lapply(shapeTerms, function(term){
    start <- sum(vapply(seq_len(term-1L), function(index)
      length(result$covPar[[index]]), integer(1)))
    start + which(result$covStruct[[term]]$free[-1L]) + 1L
  }), use.names=FALSE)
  if(length(shapeRows)) result$theta_se[shapeRows,] <- NA_real_
  if(length(shapeRows)) result$theta_se[,shapeRows] <- NA_real_
  result$convergence <- isTRUE(result$convergence) && optimum$convergence == 0L
  result
}


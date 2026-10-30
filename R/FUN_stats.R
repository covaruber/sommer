postPEV <- function(object, mode = 1L){
  
  if(!inherits(object, "mmes")){
    stop("'object' must inherit from class 'mmes'.", call. = FALSE)
  }
  
  mode <- as.integer(mode)
  
  if(length(mode) != 1L || is.na(mode) || !mode %in% 0:2){
    stop("'mode' must be one of 0, 1, or 2.", call. = FALSE)
  }

  if(inherits(object, "mmes.glmm.laplace")){
    settings <- object$laplacePostPEVSettings
    object <- .mmes_laplace_postPEV(object, mode=mode,
      panelSize=if(is.null(settings$panelSize)) 64L else settings$panelSize,
      denseMemoryMB=if(is.null(settings$denseMemoryMB)) 1000 else settings$denseMemoryMB)
    return(object)
  }

  .mmes_require_coefficient_matrix(object)

  if(!is.null(object$rotation) && mode == 1L){
    stop(
      paste0(
        "Rotated effects require mode=2 so their full eigen-coordinate ",
        "covariance can be transformed back exactly."
      ),
      call.=FALSE
    )
  }
  
  object <- post_mme_Cinverse_cpp(object, mode)

  if(!is.null(object$rotation) && mode == 2L){
    term <- object$rotation$term
    U <- object$rotation$vectors
    ranges <- object$partitions[[term]]
    for(j in seq_len(nrow(ranges))){
      rr <- ranges[j,1]:ranges[j,2]
      block <- U %*% as.matrix(object$Ci[rr,rr,drop=FALSE]) %*% t(U)
      object$uPevList[[term]][,j] <- diag(block)
    }
  }
  
  if(length(object$uPevList) && length(object$uList)){
    names(object$uPevList) <- names(object$uList)
    
    for(i in seq_along(object$uPevList)){
      if(length(object$uPevList[[i]])){
        dimnames(object$uPevList[[i]]) <- dimnames(object$uList[[i]])
      }
    }
  }
  
  object
}

.mmes_laplace_postPEV <- function(object, mode=1L, panelSize=64L,
                                  denseMemoryMB=1000){
  if(length(mode) != 1L || is.na(mode) || !mode %in% 0:2){
    stop("mode must be 0, 1, or 2.", call.=FALSE)
  }
  if(mode == 0L){
    object$Ci <- NULL
    object$VarU <- NULL
    object$uPevList <- NULL
    object$uPevMode <- 0L
    return(object)
  }
  if(!is.numeric(panelSize) || length(panelSize) != 1L ||
     !is.finite(panelSize) || panelSize < 1 || panelSize != floor(panelSize)){
    stop("output$panelSize must be a positive integer.", call.=FALSE)
  }
  if(!is.numeric(denseMemoryMB) || length(denseMemoryMB) != 1L ||
     !is.finite(denseMemoryMB) || denseMemoryMB <= 0){
    stop("execution$denseMemoryMB must be a positive finite number.", call.=FALSE)
  }
  fixedDesign <- object$fixedDesign
  randomDesign <- object$randomDesign
  curvature <- object$laplaceObservationCurvature
  randomHessian <- object$laplaceModeHessian
  fixedCount <- ncol(fixedDesign)
  randomCount <- ncol(randomDesign)
  jointCount <- fixedCount + randomCount
  curvatureDimension <- dim(curvature)
  curvatureValid <- if(is.null(curvatureDimension))
    length(curvature) == nrow(fixedDesign) else
      identical(curvatureDimension, c(nrow(fixedDesign), nrow(fixedDesign)))
  if(!curvatureValid ||
     !identical(dim(randomHessian), c(randomCount, randomCount))){
    stop("Laplace post-fit uncertainty inputs are incomplete; refit the model before requesting PEV.",
         call.=FALSE)
  }

  ordinalUncertainty <- if(.laplace_is_ordinal(object))
    .laplace_ordinal_mode_uncertainty(object, covariance=mode == 2L,
      panelSize=panelSize, denseMemoryMB=denseMemoryMB) else NULL
  if(!is.null(ordinalUncertainty)){
    randomPosteriorDiagonal <- ordinalUncertainty$diagonal
  }else{
    if(is.null(dim(curvature))){
      if(length(curvature) != nrow(fixedDesign)){
        stop("Laplace post-fit uncertainty inputs are incomplete; refit the model before requesting PEV.",
          call.=FALSE)
      }
      weightMatrix <- Matrix::Diagonal(x=curvature)
    }else{
      if(!identical(dim(curvature), c(nrow(fixedDesign), nrow(fixedDesign)))){
        stop("Laplace observation curvature has incompatible dimensions; refit the model before requesting PEV.",
          call.=FALSE)
      }
      weightMatrix <- Matrix::Matrix(curvature, sparse=TRUE)
    }
    fixedHessian <- Matrix::crossprod(fixedDesign, weightMatrix %*% fixedDesign)
    fixedRandomHessian <- Matrix::crossprod(fixedDesign, weightMatrix %*% randomDesign)
    jointHessian <- rbind(
      cbind(fixedHessian, fixedRandomHessian),
      cbind(Matrix::t(fixedRandomHessian), randomHessian))
    jointHessian <- Matrix::forceSymmetric(jointHessian)
    jointFactor <- Matrix::Cholesky(jointHessian, LDL=FALSE, perm=TRUE)
    panelSize <- min(as.integer(panelSize), jointCount)
    panelBytes <- as.double(jointCount) * panelSize * 16
    if(!is.finite(panelBytes) || panelBytes > denseMemoryMB * 1048576){
      stop("Laplace PEV solve panel exceeds execution$denseMemoryMB; reduce output$panelSize or increase the memory budget.",
           call.=FALSE)
    }
    randomPosteriorDiagonal <- numeric(randomCount)
    randomPositions <- fixedCount + seq_len(randomCount)
    for(panelStart in seq.int(1L, randomCount, by=panelSize)){
      panelStop <- min(randomCount, panelStart + panelSize - 1L)
      selected <- randomPositions[panelStart:panelStop]
      rhs <- Matrix::sparseMatrix(i=selected,
        j=seq_along(selected), x=1,
        dims=c(jointCount, length(selected)))
      solved <- Matrix::solve(jointFactor, rhs, system="A")
      randomPosteriorDiagonal[panelStart:panelStop] <-
        solved[cbind(selected, seq_along(selected))]
    }
  }

  blockWidths <- object$randomBlockWidths
  if(length(blockWidths) != length(object$covStruct) ||
     sum(blockWidths) != randomCount ||
     length(object$randomRelationships) != length(blockWidths)){
    stop("Laplace random-effect block metadata are inconsistent; refit the model before requesting PEV.",
         call.=FALSE)
  }
  termStart <- 1L
  termRanges <- vector("list", length(blockWidths))
  for(termIndex in seq_along(blockWidths)){
    relationshipCount <- nrow(object$randomRelationships[[termIndex]])
    covariance <- evaluate_covstruct_cpp(object$covStruct[[termIndex]],
      object$covStruct[[termIndex]]$par)
    coordinateCount <- nrow(covariance)
    expectedWidth <- relationshipCount * coordinateCount
    if(blockWidths[termIndex] != expectedWidth){
      stop("Laplace random-effect covariance block has an unexpected width; refit the model before requesting PEV.",
           call.=FALSE)
    }
    termStop <- termStart + blockWidths[termIndex] - 1L
    termRanges[[termIndex]] <- seq.int(termStart, termStop)
    termStart <- termStop + 1L
  }
  randomPEVDiagonal <- randomPosteriorDiagonal
  termNames <- object$randomMetadata
  if(length(termNames) != length(termRanges)){
    termNames <- paste0("random", seq_along(termRanges))
  }
  relationshipNames <- lapply(object$randomRelationships, rownames)
  if(any(vapply(relationshipNames, is.null, logical(1)))){
    relationshipNames <- lapply(object$randomRelationships, function(precision)
      as.character(seq_len(nrow(precision))))
  }
  object$uPevList <- lapply(seq_along(termRanges), function(termIndex){
    relationshipCount <- nrow(object$randomRelationships[[termIndex]])
    coordinateCount <- blockWidths[termIndex] / relationshipCount
    start <- min(termRanges[[termIndex]])
    stop <- max(termRanges[[termIndex]])
    matrix(randomPEVDiagonal[start:stop], nrow=relationshipCount,
      ncol=coordinateCount,
      dimnames=list(relationshipNames[[termIndex]],
        paste0(termNames[termIndex], ".", seq_len(coordinateCount))))
  })
  names(object$uPevList) <- termNames
  object$uPevMode <- as.integer(mode)
  object$Ci <- NULL
  object$VarU <- NULL

  if(mode == 2L && !is.null(ordinalUncertainty)){
    object$Ci <- ordinalUncertainty$covariance
  }else if(mode == 2L){
    fullBytes <- as.double(jointCount) * jointCount * 8 * 3
    if(!is.finite(fullBytes) || fullBytes > denseMemoryMB * 1048576){
      stop("Full Laplace coefficient covariance and random-effect PEV exceed execution$denseMemoryMB; use postPEV(mode=1) or increase the memory budget.",
           call.=FALSE)
    }
    jointCovariance <- as.matrix(Matrix::solve(jointFactor,
      Matrix::Diagonal(jointCount), system="A"))
    dimnames(jointCovariance) <- list(
      c(names(object$b), paste0("u", seq_len(randomCount))),
      c(names(object$b), paste0("u", seq_len(randomCount))))
    object$Ci <- jointCovariance
  }
  object
}

.mmes_laplace_prediction_uncertainty <- function(object, design,
                                                 covariance=FALSE){
  if(!isTRUE(object$convergence) ||
     length(object$laplaceDiagnostics$boundaryParameters)){
    stop("Laplace prediction uncertainty is unavailable for a nonconverged or near-boundary fit.",
      call.=FALSE)
  }
  fixedDesign <- object$fixedDesign
  randomDesign <- object$randomDesign
  curvature <- object$laplaceObservationCurvature
  randomHessian <- object$laplaceModeHessian
  fixedCount <- ncol(fixedDesign)
  randomCount <- ncol(randomDesign)
  jointCount <- fixedCount + randomCount
  curvatureDimension <- dim(curvature)
  curvatureValid <- if(is.null(curvatureDimension))
    length(curvature) == nrow(fixedDesign) else
      identical(curvatureDimension, c(nrow(fixedDesign), nrow(fixedDesign)))
  if(!curvatureValid ||
     !identical(dim(randomHessian), c(randomCount, randomCount))){
    stop("Laplace post-fit uncertainty inputs are incomplete; refit the model before requesting prediction uncertainty.",
      call.=FALSE)
  }
  weightMatrix <- if(is.null(curvatureDimension))
    Matrix::Diagonal(x=curvature) else Matrix::Matrix(curvature, sparse=TRUE)
  fixedHessian <- Matrix::crossprod(fixedDesign, weightMatrix %*% fixedDesign)
  fixedRandomHessian <- Matrix::crossprod(fixedDesign,
    weightMatrix %*% randomDesign)
  jointHessian <- rbind(
    cbind(fixedHessian, fixedRandomHessian),
    cbind(Matrix::t(fixedRandomHessian), randomHessian))
  jointHessian <- Matrix::forceSymmetric(jointHessian)
  jointFactor <- Matrix::Cholesky(jointHessian, LDL=FALSE, perm=TRUE)
  predictionCount <- nrow(design)
  if(!predictionCount){
    stop("Prediction uncertainty requires at least one prediction row.",
      call.=FALSE)
  }
  settings <- object$laplacePostPEVSettings
  denseMemoryMB <- if(is.null(settings$denseMemoryMB)) 1000 else
    settings$denseMemoryMB
  if(covariance){
    fullBytes <- as.double(jointCount) * predictionCount * 16 +
      as.double(predictionCount) * predictionCount * 8
    if(!is.finite(fullBytes) || fullBytes > denseMemoryMB * 1048576){
      stop("Full Laplace prediction covariance exceeds execution$denseMemoryMB; request se.fit=TRUE without cov.fit or increase the memory budget.",
        call.=FALSE)
    }
    solved <- Matrix::solve(jointFactor, Matrix::t(design), system="A")
    covarianceMatrix <- as.matrix(design %*% solved)
    covarianceMatrix <- (covarianceMatrix + t(covarianceMatrix)) / 2
    variance <- diag(covarianceMatrix)
    if(any(!is.finite(variance)) || any(variance < 0)){
      stop("Laplace prediction covariance has an invalid variance; no standard error was returned.",
        call.=FALSE)
    }
    return(list(se.fit=sqrt(variance),
      cov.fit=covarianceMatrix))
  }
  panelSize <- if(is.null(settings$panelSize)) 64L else settings$panelSize
  panelSize <- min(as.integer(panelSize), predictionCount)
  panelBytes <- as.double(jointCount) * panelSize * 16
  if(!is.finite(panelBytes) || panelBytes > denseMemoryMB * 1048576){
    stop("Laplace prediction uncertainty solve panel exceeds execution$denseMemoryMB; reduce output$panelSize or increase the memory budget.",
      call.=FALSE)
  }
  variance <- numeric(predictionCount)
  for(panelStart in seq.int(1L, predictionCount, by=panelSize)){
    panelStop <- min(predictionCount, panelStart + panelSize - 1L)
    rows <- panelStart:panelStop
    panelDesign <- design[rows, , drop=FALSE]
    solved <- Matrix::solve(jointFactor, Matrix::t(panelDesign), system="A")
    variance[rows] <- diag(as.matrix(panelDesign %*% solved))
  }
  if(any(!is.finite(variance)) || any(variance < 0)){
    stop("Laplace prediction variance is invalid; no standard error was returned.",
      call.=FALSE)
  }
  list(se.fit=sqrt(variance), cov.fit=NULL)
}

.mmes_laplace_population_uncertainty <- function(object, fixedDesign, offset,
                                                  randomDesign, familyIndex,
                                                  covariance=FALSE){
  uncertainty <- .mmes_laplace_uncertainty(object)
  if(!uncertainty$valid){
    stop("Population prediction uncertainty is unavailable: ",
      uncertainty$reason, ".", call.=FALSE)
  }
  parameters <- object$optimizer$par
  layout <- object$laplaceParameterLayout
  fixedCount <- layout$fixedCount
  randomWidths <- object$randomBlockWidths
  evaluateMean <- function(proposal){
    structures <- object$covStruct
    for(termIndex in seq_along(layout$freePositions)){
      freePositions <- layout$freePositions[[termIndex]]
      if(length(freePositions)){
        parameterPositions <- fixedCount + layout$parameterStarts[termIndex] +
          seq_along(freePositions)
        structures[[termIndex]]$par[freePositions] <-
          proposal[parameterPositions]
      }
    }
    linearPredictor <- offset +
      as.numeric(fixedDesign %*% proposal[seq_len(fixedCount)])
    randomVariance <- .laplace_population_variance(randomDesign, structures,
      object$randomRelationships, randomWidths)
    .laplace_population_mean(object, linearPredictor, familyIndex,
      randomVariance)
  }
  basePrediction <- evaluateMean(parameters)
  jacobian <- matrix(0, nrow=length(basePrediction),
    ncol=length(parameters))
  for(parameterIndex in seq_along(parameters)){
    step <- 1e-5 * max(1, abs(parameters[parameterIndex]))
    upper <- parameters
    lower <- parameters
    upper[parameterIndex] <- upper[parameterIndex] + step
    lower[parameterIndex] <- lower[parameterIndex] - step
    jacobian[, parameterIndex] <-
      (evaluateMean(upper) - evaluateMean(lower)) / (2 * step)
  }
  weightedJacobian <- jacobian %*% uncertainty$parameterCovariance
  variance <- rowSums(weightedJacobian * jacobian)
  if(any(!is.finite(variance)) || any(variance < 0)){
    stop("Population prediction variance is invalid; no standard error was returned.",
      call.=FALSE)
  }
  if(covariance){
    denseMemoryMB <- object$laplacePostPEVSettings$denseMemoryMB
    fullBytes <- as.double(nrow(jacobian))^2 * 8 +
      2 * as.double(nrow(jacobian)) * ncol(jacobian) * 8
    if(!is.finite(fullBytes) || fullBytes > denseMemoryMB * 1048576){
      stop("Full population prediction covariance exceeds execution$denseMemoryMB; request se.fit=TRUE without cov.fit or increase the memory budget.",
        call.=FALSE)
    }
    covarianceMatrix <- weightedJacobian %*% t(jacobian)
    covarianceMatrix <- (covarianceMatrix + t(covarianceMatrix)) / 2
    return(list(se.fit=sqrt(variance), cov.fit=covarianceMatrix))
  }
  list(se.fit=sqrt(variance), cov.fit=NULL)
}


.mmes_laplace_zero_prediction_uncertainty <- function(object, fixedDesign,
    offset, randomDesign, familyIndex, specifications, population=FALSE,
    conditional=FALSE, covariance=FALSE,
    relationships=object$randomRelationships, blockWidths=object$randomBlockWidths,
    response=TRUE, ordinalType=NULL){
  uncertainty <- .mmes_laplace_uncertainty(object)
  if(!uncertainty$valid){
    stop("Component prediction uncertainty is unavailable: ",
      uncertainty$reason, ".", call.=FALSE)
  }
  parameters <- object$optimizer$par
  layout <- object$laplaceParameterLayout
  n <- nrow(fixedDesign)
  ordinalProbabilities <- !is.null(ordinalType) && ordinalType == "response"
  outputRows <- if(ordinalProbabilities) n * length(object$ordinal$levels) else n
  budget <- object$laplacePostPEVSettings$denseMemoryMB * 1048576
  bytes <- 8 * (4 * as.double(outputRows) * length(parameters) +
    if(covariance) 3 * as.double(outputRows)^2 else 0)
  if(!is.finite(bytes) || bytes > budget){
    stop("Component prediction uncertainty exceeds execution$denseMemoryMB.",
      call.=FALSE)
  }
  context <- if(conditional) sommer_laplace_objective_create(object$fixedDesign,
    object$randomDesign, object$laplaceNativeResponse, object$offset, "",
    object$covStruct, object$randomRelationships,
    object$laplacePostPEVSettings$denseMemoryMB,
    familySpecs=object$familySpecifications, familyIndex=object$laplaceFamilyIndex,
    trials=object$priorWeights) else NULL
  evaluateMean <- function(proposal){
    local <- object
    local$familySpecifications <- .laplace_proposal_specs(specifications, proposal, layout)
    local$laplaceFamilies <- .laplace_update_families(object$laplaceFamilies,
      local$familySpecifications, familyIndex)
    eta <- offset + as.numeric(fixedDesign %*% proposal[seq_len(layout$fixedCount)])
    if(conditional){
      evaluation <- sommer_laplace_objective_eval(context, proposal,
        maxModeIterations=object$laplacePostPEVSettings$modeMaxit,
        tolerance=object$laplacePostPEVSettings$modeTol, gradient=FALSE)
      if(!isTRUE(evaluation$valid)){
        stop("Conditional prediction uncertainty encountered an invalid mode: ",
          evaluation$reason, ".", call.=FALSE)
      }
      eta <- eta + as.numeric(randomDesign %*% evaluation$mode)
    }
    if(!response || identical(ordinalType, "link")) return(eta)
    if(population){
      structures <- object$covStruct
      for(term in seq_along(layout$freePositions)){
        positions <- layout$freePositions[[term]]
        if(length(positions)) structures[[term]]$par[positions] <-
          proposal[layout$fixedCount + layout$parameterStarts[term] + seq_along(positions)]
      }
      variance <- .laplace_population_variance(randomDesign, structures,
        relationships, blockWidths)
      if(!is.null(ordinalType)) as.numeric(.laplace_ordinal_output(local, eta,
        ordinalType, randomVariance=variance)) else
        .laplace_population_mean(local, eta, familyIndex, variance)
    }else if(!is.null(ordinalType)) as.numeric(.laplace_ordinal_output(local, eta,
      ordinalType)) else .pql_apply(local$laplaceFamilies, familyIndex, "linkinv", eta)
  }
  jacobian <- matrix(0, outputRows, length(parameters))
  for(position in seq_along(parameters)){
    step <- 1e-5 * max(1, abs(parameters[position]))
    upper <- lower <- parameters
    upper[position] <- upper[position] + step
    lower[position] <- lower[position] - step
    jacobian[, position] <- (evaluateMean(upper) - evaluateMean(lower)) / (2 * step)
  }
  weighted <- jacobian %*% uncertainty$parameterCovariance
  variance <- rowSums(weighted * jacobian)
  resultCovariance <- if(covariance) weighted %*% t(jacobian) else NULL
  if(conditional){
    if(16 * as.double(ncol(randomDesign)) * outputRows > budget){
      stop("Conditional zero-component uncertainty solve exceeds execution$denseMemoryMB.",
        call.=FALSE)
    }
    factor <- Matrix::Cholesky(Matrix::forceSymmetric(object$laplaceModeHessian),
      LDL=FALSE)
    solved <- Matrix::solve(factor, Matrix::t(randomDesign), system="A")
    eta <- offset + as.numeric(fixedDesign %*% object$b) +
      as.numeric(randomDesign %*% object$randomMode)
    families <- .laplace_update_families(object$laplaceFamilies, specifications, familyIndex)
    derivative <- if(!is.null(ordinalType)){
      if(ordinalType == "link") rep(1, n) else {
        distribution <- .laplace_ordinal_distribution(eta, specifications[[1L]],
          object$ordinal$levels, object$laplacePostPEVSettings$denseMemoryMB)
        if(ordinalProbabilities) as.numeric(distribution$derivative) else
          as.numeric(distribution$derivative %*% seq_len(ncol(distribution$derivative)))
      }
    }else if(response) .pql_apply(families, familyIndex, "mu.eta", eta) else rep(1, n)
    if(ordinalProbabilities){
      rows <- rep(seq_len(n), length(object$ordinal$levels))
      randomDesign <- randomDesign[rows, , drop=FALSE]
      solved <- solved[, rows, drop=FALSE]
    }
    variance <- variance + derivative^2 *
      as.numeric(Matrix::rowSums(randomDesign * Matrix::t(solved)))
    if(covariance) resultCovariance <- resultCovariance +
      outer(derivative, derivative) * as.matrix(randomDesign %*% solved)
  }
  if(any(!is.finite(variance)) || any(variance < 0)){
    stop("Zero-component prediction variance is invalid; no standard error was returned.",
      call.=FALSE)
  }
  list(se.fit=sqrt(variance), cov.fit=resultCovariance)
}

.mmes_varu_check <- function(object){
  if(!inherits(object, "mmes") || inherits(object, "mmes.glmm")){
    stop("VarU requires a Gaussian mmes fit.", call.=FALSE)
  }
  .mmes_require_coefficient_matrix(object)
  if(is.null(object$C) || !nrow(object$C) || is.null(object$Cscale)){
    stop("VarU requires a Henderson fit with stored C and Cscale; direct and task$mode='solve' fits are not supported.", call.=FALSE)
  }
}

.mmes_random_precision <- function(object){
  precision <- object$randomPrecision
  if(is.null(precision)){
    inputs <- .mmes_rebuild_inputs(object)
    if(!is.null(object$rotation) || !is.null(object$factorScoreInfo)){
      stop("Refit rotated or factor-score models to retain original relationship precisions for VarU.", call.=FALSE)
    }
    precision <- inputs$Ai[seq_along(object$uList)]
    for(index in seq_along(precision)){
      if(nrow(precision[[index]]) != nrow(object$uList[[index]])){
        stop("Rebuilt relationship precision does not match fitted effects; refit the model.", call.=FALSE)
      }
    }
  }
  precision
}

.mmes_prior_vcov <- function(object, D, precision=.mmes_random_precision(object),
                             factors=NULL){
  out <- matrix(0, nrow(D), nrow(D))
  for(index in seq_along(object$uList)){
    ranges <- object$partitions[[index]]
    columns <- unlist(lapply(seq_len(nrow(ranges)), function(coordinate){
      seq.int(ranges[coordinate,1], ranges[coordinate,2])
    }), use.names=FALSE)
    if(!any(D[,columns,drop=FALSE] != 0)) next
    Sigma <- object$theta[[index]]
    dimension <- object$covStruct[[index]]$dim
    if(!identical(dim(Sigma), c(dimension, dimension))){
      Sigma <- covmatrix_mmes(object, index, se=FALSE, max.dim=dimension)$covariance
    }
    blocks <- lapply(seq_len(nrow(ranges)), function(coordinate){
      D[,seq.int(ranges[coordinate,1], ranges[coordinate,2]),drop=FALSE]
    })
    factor <- if(is.null(factors)){
      Matrix::Cholesky(Matrix::forceSymmetric(precision[[index]]), LDL=FALSE)
    }else factors[[index]]
    solved <- Matrix::solve(factor, do.call(cbind, lapply(blocks, t)), system="A")
    for(first in seq_along(blocks)){
      for(second in seq_along(blocks)){
        rhsColumns <- (second-1L)*nrow(D) + seq_len(nrow(D))
        out <- out + Sigma[first,second] * as.matrix(blocks[[first]] %*% solved[,rhsColumns,drop=FALSE])
      }
    }
  }
  (out+t(out))/2
}

.mmes_vcov_solver <- function(object){
  .mmes_require_coefficient_matrix(object)
  factor <- mmes_vcov_factor_cpp(object)
  function(D){
    mmes_vcov_solve_cpp(factor, to_sparse(.mmes_engine_contrast(object, D)))
  }
}

.mmes_henderson_vcov <- function(object, contrasts){
  solve <- .mmes_vcov_solver(object)
  lapply(contrasts, solve)
}

postVarU <- function(object, mode=1L){
  .mmes_varu_check(object)
  if(length(mode) != 1L || is.na(mode) || !mode %in% 0:2){
    stop("mode must be 0, 1, or 2.", call.=FALSE)
  }
  if(mode == 0L){
    object$uVarList <- NULL
    object$VarU <- NULL
    object$uVarMode <- 0L
    return(object)
  }
  precision <- .mmes_random_precision(object)
  priorFactors <- lapply(precision, function(x)
    Matrix::Cholesky(Matrix::forceSymmetric(x), LDL=FALSE))
  solve <- .mmes_vcov_solver(object)
  nFixed <- nrow(object$b)
  nRandom <- nrow(object$bu)-nFixed
  diagonal <- numeric(nRandom)
  if(mode == 2L){
    D <- cbind(Matrix::Matrix(0, nRandom, nFixed, sparse=TRUE), Matrix::Diagonal(nRandom))
    error <- solve(D)
    object$VarU <- .mmes_prior_vcov(object, D, precision, priorFactors)-error
    labels <- unlist(lapply(seq_along(object$uList), function(index){
      effects <- object$uList[[index]]
      paste(names(object$uList)[index],
            rep(colnames(effects), each=nrow(effects)),
            rep(rownames(effects), times=ncol(effects)), sep=":")
    }))
    dimnames(object$VarU) <- list(labels, labels)
    diagonal <- diag(object$VarU)
  }else{
    object$VarU <- NULL
    for(start in seq.int(1L, max(nRandom,1L), by=128L)){
      if(!nRandom) break
      rows <- seq.int(start, min(start+127L, nRandom))
      D <- Matrix::sparseMatrix(i=seq_along(rows), j=nFixed+rows, x=1,
                                dims=c(length(rows), nrow(object$bu)))
      error <- solve(D)
      diagonal[rows] <- diag(.mmes_prior_vcov(object, D, precision, priorFactors)-error)
    }
  }
  object$uVarList <- lapply(seq_along(object$uList), function(index){
    ranges <- object$partitions[[index]]
    columns <- unlist(lapply(seq_len(nrow(ranges)), function(coordinate){
      seq.int(ranges[coordinate,1], ranges[coordinate,2])
    }), use.names=FALSE)-nFixed
    matrix(diagonal[columns], nrow=nrow(object$uList[[index]]),
           dimnames=dimnames(object$uList[[index]]))
  })
  names(object$uVarList) <- names(object$uList)
  object$uVarMode <- as.integer(mode)
  object
}

corImputation <- function(wide, Gu=NULL, nearest=10, roundR=FALSE){
  if(is.null(rownames(wide))){stop("Rownames of the input matrix cannot be NULL. Please add them", call. = FALSE)}
  if(is.null(Gu)){
    X <- apply(wide, 2, imputev)
    Gu <- cor(t(X))
  }
  wide2 <- wide
  rowNamesWide <-  rownames(wide)
  # for each feature
  for(iEnv in 1:ncol(wide)){ # iEnv=10
    withData <- which(!is.na(wide[,iEnv]))
    withoutData <- which(is.na(wide[,iEnv]))
    toPredict <- 1:nrow(wide)
    # if(length(toPredict) > 1){
    likelihood=Gu[as.character(rowNamesWide)[toPredict],as.character(rowNamesWide)[withData], drop=FALSE]
    # }else{
    #   likelihood=matrix(Gu[as.character(rowNamesWide)[toPredict],as.character(rowNamesWide)[withData]], nrow=1)
    #   rownames(likelihood) <- as.character(rowNamesWide)[toPredict]
    #   colnames(likelihood) <- as.character(rowNamesWide)[withData]
    # }
    replacement <- numeric()
    for(iInd in 1:nrow(likelihood)){ # iInd=1
      # averaging only the positively correlated to avoid decrease   # wide[iInd, iEnv]
      indLik <- sort(abs(likelihood[iInd,]), decreasing = TRUE)
      toAverage <- indLik[1:min(c(nearest,length(withData)))]
      indLikToAverage <- likelihood[iInd,names(toAverage)]
      replacement[iInd] <- mean(wide[names(which(indLikToAverage > 0)),iEnv]) 
    }
    names(replacement) <- rownames(likelihood)
    # time to replace the missing data
    dd=data.frame(replacement=replacement, index=1:length(replacement),
                  imputed=1, id=rownames(likelihood),orVal=wide[,iEnv])
    dd$imputed[which(dd$id %in% names(withoutData))]=0
    head(dd)
    if(roundR){
      wide2[toPredict,iEnv] <- round(replacement)
    }else{
      wide2[toPredict,iEnv] <- replacement
    }
  }
  # for each individual
  for(jRow in 1:nrow(wide)){
    miss <- which(is.na(wide[jRow,]))
    if(length(miss) > 0){
      dd=data.frame(full=as.vector(unlist(wide2[jRow,])), partial=as.vector(unlist(wide[jRow,])))
      # model <- fastLm(partial~full,data=dd[which(!is.na(dd$partial)),])
      # model <- lm(partial~full,data=dd[which(!is.na(dd$partial)),])
      model <- mmes(partial~full,data=dd[which(!is.na(dd$partial)),],
                    control=mmes_control(diagnostics=list(verbose=FALSE)))
      pp=as.vector(model$b[1,1])+(dd[which(is.na(dd$partial)),"full"]*as.vector(model$b[2,1]))
      if(roundR){
        wide[jRow,miss] <- round(pp)#round(predict(model,newdata = dd[which(is.na(dd$partial)),]))
      }else{
        wide[jRow,miss] <- pp#predict(model,newdata = dd[which(is.na(dd$partial)),])
      }
    }
  }
  stillEmpty <- which(is.na(wide), arr.ind = TRUE)
  if(nrow(stillEmpty) > 0){wide[stillEmpty] <- mean(wide, na.rm=TRUE)}
  stillEmpty <- which(is.na(wide2), arr.ind = TRUE)
  if(nrow(stillEmpty) > 0){wide2[stillEmpty] <- mean(wide2, na.rm=TRUE)}
  
  return(list(imputed=wide, corImputed=wide2))
}


r2 <- function(object, object2=NULL){
  if(!inherits(object, "mmes")){
    stop("This function is only available for models fitted with the mmes() function.", call. = FALSE)
  }
  result <- list()
  for(iPart in 1:length(object$uPevList)){
    pev <- object$uPevList[[iPart]]
    variance0 <- as.numeric(diag(object$theta[[iPart]]))
    variance <- apply(data.frame(variance0),1,function(x){rep(x,nrow(pev))})
    if(!is.null(object2)){
      if(!is.null(object2$Ai)){
        relationshipDiagonal <- diag(solve(object2$Ai[[iPart]]))
        for(iVar in 1:length(variance0)){ # if user provided A matrices we subsitute with more accurate values
          variance[,iVar] <- variance0[iVar] * relationshipDiagonal
        }
      }else{
        stop("object2 needs to be a model that sets the argument 'returnParam=TRUE' so we can extract the relationship matrices. Please correct. ", call. = FALSE)
      }
      
    }
    result[[iPart]] <- (variance - pev)/variance
  }
  names(result) <- names(object$uList)
  return(result)
}




matrix.trace <- function(x){
  if (!is.square.matrix(x))
    stop("argument x is not a square matrix")
  return(sum(diag(x)))
}

is.diagonal.matrix <- function (x, tol = 1e-08){
  y <- x
  diag(y) <- rep(0, nrow(y))
  return(all(abs(y) < tol))
}

is.square.matrix <-function(x){
  return(nrow(x) == ncol(x))
}

hadamard.prod <-function (x, y){
  if (!is.numeric(x)) {
    stop("argument x is not numeric")
  }
  if (!is.numeric(y)) {
    stop("argument y is not numeric")
  }
  if (is.matrix(x)) {
    Xmat <- x
  }
  else {
    if (is.vector(x)) {
      Xmat <- matrix(x, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (is.matrix(y)) {
    Ymat <- y
  }
  else {
    if (is.vector(y)) {
      Ymat <- matrix(y, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (nrow(Xmat) != nrow(Ymat))
    stop("argumentx x and y do not have the same row order")
  if (ncol(Xmat) != ncol(Ymat))
    stop("arguments x and y do not have the same column order")
  return(Xmat * Ymat)
}





# Least squares for variance-component starting values only. A well-conditioned
# full-rank design goes through the normal equations (sparse crossprod plus a
# Cholesky factor); anything else uses the dense QR, whose NA coefficients flag
# aliased columns to the caller.
.mmes_start_lsq <- function(X, Y){
  p <- ncol(X)
  XtX <- as.matrix(Matrix::crossprod(X))
  R <- tryCatch(chol(XtX), error=function(e) NULL)
  if(!is.null(R) && all(is.finite(R)) &&
     min(diag(R) / sqrt(diag(XtX))) > 1e-7){
    XtY <- as.matrix(Matrix::crossprod(X, Y))
    beta <- backsolve(R, backsolve(R, XtY, transpose=TRUE))
    return(list(beta=beta, rank=p, resid=Y - as.matrix(X %*% beta)))
  }
  Xd <- as.matrix(X)
  qrX <- qr(Xd)
  beta <- if(nrow(Xd) - qrX$rank >= 1L) qr.coef(qrX, Y) else NULL
  list(beta=beta, rank=qrX$rank,
       resid=if(is.null(beta) || anyNA(beta)) NULL else Y - Xd %*% beta)
}

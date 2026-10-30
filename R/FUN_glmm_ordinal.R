.laplace_ordinal_configuration <- function(family, settings){
  if(is.null(settings)) settings <- list()
  configuration <- .mme_module(settings,
    list(thresholds=family$thresholds, estimate=TRUE),
    "familyParameters")
  configuration$estimate <- .mme_flag(configuration$estimate,
    "familyParameters$estimate")
  # Apply the same cutpoint validation to constructor and control inputs.
  configuration$thresholds <- ordinal(family$link,
    configuration$thresholds)$thresholds
  if(!configuration$estimate && is.null(configuration$thresholds)){
    stop("Fixed ordinal cutpoints require familyParameters$thresholds or ordinal(thresholds=).",
      call.=FALSE)
  }
  configuration
}

.laplace_ordinal_coordinates <- function(thresholds){
  thresholds <- ordinal(thresholds=thresholds)$thresholds
  if(is.null(thresholds)){
    stop("Ordinal cutpoint coordinates require a nonempty threshold vector.",
      call.=FALSE)
  }
  gaps <- diff(thresholds)
  if(any(!is.finite(gaps))){
    stop("Ordinal cutpoint gaps exceed the finite numerical range.", call.=FALSE)
  }
  c(thresholds[1L], log(gaps))
}

.laplace_ordinal_thresholds <- function(coordinates){
  if(!is.numeric(coordinates) || is.complex(coordinates) || !is.null(dim(coordinates)) ||
     !length(coordinates) || any(!is.finite(coordinates))){
    stop("Ordinal cutpoint coordinates must be a nonempty finite numeric vector.",
      call.=FALSE)
  }
  thresholds <- cumsum(c(coordinates[1L], exp(coordinates[-1L])))
  if(any(!is.finite(thresholds)) || any(diff(thresholds) <= 0)){
    stop("Ordinal cutpoint proposal is outside the ordered finite numerical range.",
      call.=FALSE)
  }
  thresholds
}

.laplace_prepare_ordinal <- function(fixed, random, rcov, data, family,
                                     naMethodY, mmesArgs, familyParameters){
  configuration <- .laplace_ordinal_configuration(family, familyParameters)
  frame <- stats::model.frame(fixed, data=data, na.action=stats::na.pass,
    drop.unused.levels=FALSE)
  response <- stats::model.response(frame)
  if(!is.null(dim(frame[[1L]])) ||
     (!is.ordered(response) && (!is.numeric(response) ||
       is.object(response) || is.complex(response)))){
    stop("Ordinal responses must be ordered factors or plain numeric integer category codes; unordered factors and characters are not supported.",
      call.=FALSE)
  }
  if(!is.null(attr(response, "trials"))){
    stop("Ordinal responses do not support binomial trials.", call.=FALSE)
  }
  orderedResponse <- is.ordered(response)
  codes <- if(orderedResponse) as.integer(response) else response
  observed <- codes[!is.na(codes)]
  if(!length(observed) || any(!is.finite(observed)) ||
     any(observed < 1 | observed != floor(observed))){
    stop("Ordinal category codes must be finite positive integers with observed categories.",
      call.=FALSE)
  }
  labels <- if(orderedResponse) levels(response) else
    as.character(sort(unique(observed)))
  categoryCount <- length(labels)
  if(categoryCount < 2L || (!orderedResponse &&
     !identical(sort(unique(as.numeric(observed))),
       as.numeric(seq_len(categoryCount))))){
    stop("Ordinal responses require at least two categories; numeric category codes must be contiguous starting at one.",
      call.=FALSE)
  }
  thresholds <- configuration$thresholds
  if(!is.null(thresholds) && length(thresholds) != categoryCount - 1L){
    stop("Ordinal responses require exactly one fewer cutpoint than declared categories.",
      call.=FALSE)
  }

  numericFixed <- fixed
  if(orderedResponse){
    numericFixed[[2L]] <- substitute(base::as.integer(responseExpression),
      list(responseExpression=fixed[[2L]]))
  }
  prepared <- .mmes_prepare_model(numericFixed, random, rcov, data, family,
    naMethodY, mmesArgs, .ordinalPrepared=TRUE)
  included <- prepared$obsInfo$included
  counts <- tabulate(codes[included], nbins=categoryCount)
  if(!sum(counts)){
    stop("Ordinal preparation requires retained observations.", call.=FALSE)
  }
  if(configuration$estimate && any(counts == 0L)){
    stop("Estimated ordinal cutpoints require retained observations in every declared category.",
      call.=FALSE)
  }
  if(is.null(thresholds)){
    cumulative <- head(cumsum(counts) / sum(counts), -1L)
    thresholds <- family$linkfun(cumulative)
  }
  coordinates <- .laplace_ordinal_coordinates(thresholds)
  design <- prepared$X
  fullColumns <- colnames(design)
  fittedContrasts <- attr(design, "contrasts")
  assignment <- attr(design, "assign")
  if(is.null(assignment)){
    assignment <- attr(stats::model.matrix(attr(prepared$modelFrame, "terms"),
      prepared$retainedFrame, contrasts.arg=mmesArgs$contrasts), "assign")
  }
  columns <- seq_len(ncol(design))
  if(configuration$estimate) columns <- columns[assignment != 0L]
  design <- design[, columns, drop=FALSE]
  if(ncol(design) && Matrix::rankMatrix(design)[1L] != ncol(design)){
    stop("Ordinal preparation requires a full-rank fixed-effects design.",
      call.=FALSE)
  }
  if(configuration$estimate &&
     Matrix::rankMatrix(cbind(design, rep(1, nrow(design))))[1L] != ncol(design) + 1L){
    stop("Estimated ordinal cutpoint locations are confounded with the fixed design after intercept removal.",
      call.=FALSE)
  }
  prepared$X <- design
  prepared$partitionsX <- lapply(prepared$partitionsX, function(partition){
    indices <- match(as.vector(partition), columns)
    matrix(indices[!is.na(indices)], nrow=1L)
  })
  prepared$partitionsX <- Filter(function(partition) length(partition) > 0L,
    prepared$partitionsX)
  prepared$ordinal <- list(levels=labels, counts=counts, response=as.numeric(codes[included]),
    originalFixed=fixed, originalResponse=response, columns=columns,
    fullColumns=fullColumns, contrasts=fittedContrasts,
    fixedTerms=stats::delete.response(attr(prepared$modelFrame, "terms")),
    xlevels=stats::.getXlevels(stats::delete.response(attr(prepared$modelFrame, "terms")),
      prepared$retainedFrame),
    assignment=assignment[columns],
    specification=list(family="ordinal", link=family$link, thresholds=thresholds,
      logParameter=0, free=configuration$estimate, parameter="thresholds",
      supplied=!is.null(configuration$thresholds)),
    initialParameters=if(configuration$estimate) coordinates else numeric())
  prepared
}

.laplace_is_ordinal <- function(object){
  inherits(object$family, "family") && identical(object$family$family, "ordinal")
}

.laplace_ordinal_distribution <- function(eta, specification, labels,
                                          maxMemoryMB=1000){
  n <- length(eta)
  k <- length(specification$thresholds) + 1L
  if(8 * as.double(n) * (4 * k + 12) > maxMemoryMB * 1048576){
    stop("Ordinal category probabilities exceed execution$denseMemoryMB.", call.=FALSE)
  }
  if(any(!is.finite(eta))){
    stop("Ordinal prediction requires finite predictors.", call.=FALSE)
  }
  probability <- derivative <- logProbability <- matrix(0, n, k,
    dimnames=list(names(eta), labels))
  specification$free <- FALSE
  for(category in seq_len(k)){
    value <- sommer_glmm_family_eval(rep(category, n), eta, list(specification),
      rep(1L, n), rep(1, n))
    logProbability[, category] <- as.numeric(value$logDensityRows)
    probability[, category] <- exp(logProbability[, category])
    derivative[, category] <- probability[, category] * as.numeric(value$score)
  }
  if(any(!is.finite(probability)) || any(!is.finite(derivative)) ||
     any(abs(rowSums(probability) - 1) > 1e-10)){
    stop("Ordinal category distribution is numerically invalid.", call.=FALSE)
  }
  list(probability=probability, derivative=derivative, logProbability=logProbability)
}

.laplace_ordinal_table <- function(specification, labels){
  data.frame(threshold=paste(head(labels, -1L), tail(labels, -1L), sep="|"),
    estimate=specification$thresholds, StdError=NA_real_,
    estimated=specification$free, stringsAsFactors=FALSE)
}

.laplace_ordinal_jacobian <- function(specification){
  n <- length(specification$thresholds)
  result <- matrix(0, n, n)
  result[, 1L] <- 1
  if(n > 1L) for(position in 2:n){
    result[position:n, position] <- diff(specification$thresholds)[position - 1L]
  }
  result
}

.laplace_ordinal_output <- function(object, eta, type, specification=object$familySpecifications[[1L]],
                                    randomVariance=NULL){
  if(type == "link") return(eta)
  distribution <- if(is.null(randomVariance)){
    .laplace_ordinal_distribution(eta, specification, object$ordinal$levels,
      object$laplacePostPEVSettings$denseMemoryMB)$probability
  }else .laplace_ordinal_population(object, eta, specification, randomVariance)
  if(type == "mean") return(as.numeric(distribution %*% seq_len(ncol(distribution))))
  distribution
}

.laplace_ordinal_population <- function(object, eta, specification, variance){
  if(8 * as.double(length(eta)) * (length(object$ordinal$levels) * 4 + 12) >
     object$laplacePostPEVSettings$denseMemoryMB * 1048576){
    stop("Ordinal population probabilities exceed execution$denseMemoryMB.", call.=FALSE)
  }
  if(specification$link == "probit"){
    # Convolution of two normal variables is exact on this link.
    result <- matrix(0, length(eta), length(specification$thresholds) + 1L,
      dimnames=list(NULL, object$ordinal$levels))
    for(row in seq_along(eta)){
      scale <- sqrt(1 + variance[row])
      local <- specification
      local$thresholds <- local$thresholds / scale
      result[row, ] <- .laplace_ordinal_distribution(eta[row] / scale, local,
        object$ordinal$levels, object$laplacePostPEVSettings$denseMemoryMB)$probability
    }

    return(result)
  }
  quadrature <- .laplace_normal_quadrature()
  result <- matrix(0, length(eta), length(object$ordinal$levels),
    dimnames=list(NULL, object$ordinal$levels))
  for(position in seq_along(quadrature$nodes)){
    result <- result + quadrature$weights[position] *
      .laplace_ordinal_distribution(eta + sqrt(variance) * quadrature$nodes[position],
        specification, object$ordinal$levels,
        object$laplacePostPEVSettings$denseMemoryMB)$probability
  }
  result
}

.laplace_ordinal_mode_uncertainty <- function(object, covariance=FALSE,
    panelSize=object$laplacePostPEVSettings$panelSize,
    denseMemoryMB=object$laplacePostPEVSettings$denseMemoryMB){
  object$laplacePostPEVSettings$panelSize <- panelSize
  object$laplacePostPEVSettings$denseMemoryMB <- denseMemoryMB
  uncertainty <- .mmes_laplace_uncertainty(object)
  if(!uncertainty$valid){
    stop("Ordinal PEV is unavailable: ", uncertainty$reason, ".", call.=FALSE)
  }
  parameters <- object$optimizer$par
  q <- length(object$randomMode)
  bytes <- 8 * (4 * as.double(q) * length(parameters) +
    if(covariance) 3 * as.double(q + length(object$b))^2 else 0)
  if(bytes > denseMemoryMB * 1048576){
    stop("Ordinal PEV exceeds execution$denseMemoryMB.", call.=FALSE)
  }
  context <- sommer_laplace_objective_create(object$fixedDesign,
    object$randomDesign, object$laplaceNativeResponse, object$offset, "",
    object$covStruct, object$randomRelationships, denseMemoryMB,
    familySpecs=object$familySpecifications, familyIndex=object$laplaceFamilyIndex,
    trials=object$priorWeights)
  evaluate <- function(parameters){
    value <- sommer_laplace_objective_eval(context, parameters, gradient=FALSE,
      maxModeIterations=object$laplacePostPEVSettings$modeMaxit,
      tolerance=object$laplacePostPEVSettings$modeTol)
    if(!value$valid) stop("Ordinal PEV encountered an invalid mode: ", value$reason,
      call.=FALSE)
    as.numeric(value$mode)
  }
  jacobian <- matrix(0, q, length(parameters))
  for(position in seq_along(parameters)){
    upper <- lower <- parameters
    step <- 1e-5 * max(1, abs(parameters[position]))
    upper[position] <- upper[position] + step
    lower[position] <- lower[position] - step
    jacobian[, position] <- (evaluate(upper) - evaluate(lower)) / (2 * step)
  }
  weighted <- jacobian %*% uncertainty$parameterCovariance
  factor <- Matrix::Cholesky(Matrix::forceSymmetric(object$laplaceModeHessian), LDL=FALSE)
  panelSize <- min(panelSize, q)
  if(16 * as.double(q) * panelSize > denseMemoryMB * 1048576){
    stop("Ordinal PEV solve panel exceeds execution$denseMemoryMB.", call.=FALSE)
  }
  diagonal <- rowSums(weighted * jacobian)
  for(first in seq.int(1L, q, by=panelSize)){
    columns <- seq.int(first, min(q, first + panelSize - 1L))
    rhs <- Matrix::sparseMatrix(i=columns, j=seq_along(columns), x=1,
      dims=c(q, length(columns)))
    solved <- Matrix::solve(factor, rhs, system="A")
    diagonal[columns] <- diagonal[columns] + solved[cbind(columns, seq_along(columns))]
  }
  joint <- NULL
  if(covariance){
    p <- length(object$b)
    map <- rbind(cbind(diag(p), matrix(0, p, length(parameters) - p)), jacobian)
    joint <- map %*% uncertainty$parameterCovariance %*% t(map)
    random <- p + seq_len(q)
    joint[random, random] <- joint[random, random] +
      as.matrix(Matrix::solve(factor, Matrix::Diagonal(q), system="A"))
    dimnames(joint) <- list(c(names(object$b), paste0("u", seq_len(q))),
      c(names(object$b), paste0("u", seq_len(q))))
  }
  if(any(!is.finite(diagonal)) || any(diagonal < 0)){
    stop("Ordinal PEV has an invalid variance; no uncertainty was returned.", call.=FALSE)
  }
  list(diagonal=diagonal, covariance=joint)
}

.laplace_ordinal_residuals <- function(object, type){
  distribution <- .laplace_ordinal_distribution(object$linear.predictors,
    object$familySpecifications[[1L]], object$ordinal$levels,
    object$laplacePostPEVSettings$denseMemoryMB)
  p <- distribution$probability
  response <- object$y
  mean <- object$fitted.values
  if(type == "response") return(response - mean)
  if(type == "quantile"){
    lower <- upper <- numeric(length(response))
    for(category in seq_len(ncol(p))){
      upper <- upper + p[, category] * (response >= category)
      lower <- lower + p[, category] * (response > category)
    }
    return(stats::qnorm(pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps,
      lower + stats::runif(length(response)) * (upper - lower)))))
  }
  if(type == "deviance"){
    return(ifelse(response >= mean, 1, -1) * sqrt(pmax(0,
      -2 * distribution$logProbability[cbind(seq_along(response), response)])))
  }
  denominator <- if(type == "pearson"){
    sqrt(pmax(0, as.numeric(p %*% seq_len(ncol(p))^2) - mean^2))
  }else as.numeric(distribution$derivative %*% seq_len(ncol(p)))
  if(any(!is.finite(denominator)) || any(denominator <= 0)){
    stop("Ordinal residuals are undefined at zero variance or derivative.", call.=FALSE)
  }
  (response - mean) / denominator
}

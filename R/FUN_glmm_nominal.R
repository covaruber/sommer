.laplace_nominal_response <- function(response, reference=NULL){
  if(!is.null(dim(response)) || is.ordered(response) ||
     (!is.factor(response) && (!is.numeric(response) || is.object(response) ||
       is.complex(response)))){
    stop("Nominal responses must be unordered factors or plain numeric integer category codes; convert ordered outcomes explicitly if their order is to be ignored.",
      call.=FALSE)
  }
  if(!is.null(attr(response, "trials"))){
    stop("Nominal responses do not support binomial trials.", call.=FALSE)
  }
  observed <- response[!is.na(response)]
  if(!length(observed)){
    stop("Nominal responses require observed categories.", call.=FALSE)
  }
  if(is.factor(response)){
    labels <- levels(response)
    codes <- as.integer(response)
  }else{
    observed <- sort(unique(observed))
    if(any(!is.finite(observed)) || any(observed < 1 | observed != floor(observed)) ||
       !identical(as.numeric(observed), as.numeric(seq_along(observed)))){
      stop("Nominal numeric category codes must be finite contiguous integers starting at one.",
        call.=FALSE)
    }
    labels <- as.character(observed)
    codes <- as.integer(response)
  }
  if(length(labels) < 2L || anyNA(labels) || any(!nzchar(labels)) ||
     anyDuplicated(labels)){
    stop("Nominal responses require at least two distinct nonempty category labels.",
      call.=FALSE)
  }
  if(is.null(reference)) reference <- tail(labels, 1L)
  if(!is.character(reference) || length(reference) != 1L || is.na(reference) ||
     !reference %in% labels){
    stop("Nominal reference must be one declared category label.", call.=FALSE)
  }
  nativeLevels <- c(labels[labels != reference], reference)
  list(levels=labels, reference=reference, nativeLevels=nativeLevels,
    predictorLevels=head(nativeLevels, -1L),
    response=as.numeric(match(labels[codes], nativeLevels)),
    originalResponse=response)
}

multinomial <- function(reference=NULL){
  if(!is.null(reference) && (!is.character(reference) || length(reference) != 1L ||
     is.na(reference) || !nzchar(reference))){
    stop("multinomial(reference=) must be NULL or one nonempty category label.",
      call.=FALSE)
  }
  structure(list(family="multinomial", link="softmax", reference=reference),
    class=c("sommer_nominal_family", "family"))
}

.laplace_is_nominal <- function(object){
  inherits(object$family, "sommer_nominal_family") &&
    !is.null(object$nominal)
}

.laplace_nominal_probabilities <- function(linearPredictors, nominal){
  if(!is.matrix(linearPredictors) ||
     ncol(linearPredictors) != length(nominal$predictorLevels) ||
     any(!is.finite(linearPredictors))){
    stop("Nominal linear predictors must be a finite observation-by-category matrix.",
      call.=FALSE)
  }
  shift <- pmax(0, apply(linearPredictors, 1L, max))
  nonreference <- exp(linearPredictors - shift)
  reference <- exp(-shift)
  denominator <- rowSums(nonreference) + reference
  nativeProbability <- cbind(nonreference, reference) / denominator
  probability <- nativeProbability[
    , match(nominal$levels, nominal$nativeLevels), drop=FALSE]
  colnames(probability) <- nominal$levels
  rownames(probability) <- rownames(linearPredictors)
  probability
}

.laplace_nominal_deviance <- function(linearPredictors, nominal, response){
  if(!is.matrix(linearPredictors) ||
     ncol(linearPredictors) != length(nominal$predictorLevels) ||
     any(!is.finite(linearPredictors)) ||
     length(response) != nrow(linearPredictors) || anyNA(response) ||
     any(!response %in% seq_along(nominal$nativeLevels))){
    stop("Nominal deviance requires finite logits and one valid category per observation.",
      call.=FALSE)
  }
  logits <- cbind(linearPredictors, 0)
  shift <- apply(logits, 1L, max)
  centered <- logits - shift
  contribution <- 2 * (log(rowSums(exp(centered))) -
    centered[cbind(seq_along(response), response)])
  if(any(!is.finite(contribution)) || any(contribution < 0)){
    stop("Nominal deviance contributions exceed the finite non-negative numerical range.",
      call.=FALSE)
  }
  names(contribution) <- rownames(linearPredictors)
  contribution
}

.laplace_nominal_residuals <- function(object, type){
  if(!type %in% c("response", "pearson", "deviance")){
    stop("Nominal Laplace residuals support response, pearson and unsigned deviance; working and quantile residuals are unsupported for unordered categories.",
      call.=FALSE)
  }
  if(type == "deviance"){
    response <- match(as.character(object$y), object$nominal$nativeLevels)
    return(sqrt(.laplace_nominal_deviance(object$linear.predictors,
      object$nominal, response)))
  }
  observed <- stats::model.matrix(~object$y - 1L)
  colnames(observed) <- levels(object$y)
  residual <- observed - object$category.probabilities
  if(type == "response") return(residual)
  probability <- object$category.probabilities
  if(any(!is.finite(probability)) || any(probability <= 0 | probability >= 1)){
    stop("Nominal Pearson residuals are undefined for numerically zero or one category probabilities.",
      call.=FALSE)
  }
  residual / sqrt(probability * (1 - probability))
}

.laplace_nominal_population_prediction <- function(object, prediction){
  d <- length(object$nominal$predictorLevels)
  n <- prediction$n
  if(d > 3L){
    stop("Nominal population prediction supports at most three non-reference logits.",
      call.=FALSE)
  }
  covariance <- .laplace_latent_covariance(prediction$random,object$covStruct,
    object$randomRelationships,object$randomBlockWidths,
    matrix(seq_len(n*d),nrow=n),
    object$laplacePostPEVSettings$denseMemoryMB)
  eta <- matrix(as.numeric(prediction$fixed %*% object$b)+prediction$offset,
    nrow=n,ncol=d)
  result <- matrix(0,n,d+1L)
  remaining <- seq_len(n)
  while(length(remaining)){
    observation <- remaining[1L]
    same <- vapply(remaining,function(index)
      identical(covariance[,,index],covariance[,,observation]),logical(1))
    selected <- remaining[same]
    result[selected,] <- .laplace_nominal_integrate(eta[selected,,drop=FALSE],
      matrix(covariance[,,observation],d,d),object$nominal,
      object$laplacePostPEVSettings$denseMemoryMB)
    remaining <- remaining[!same]
  }
  dimnames(result) <- list(prediction$rowNames,object$nominal$levels)
  result
}

.laplace_nominal_integrate <- function(eta, covariance, specification, memoryMB){
  d <- ncol(eta)
  decomposition <- eigen(covariance,symmetric=TRUE)
  if(any(!is.finite(decomposition$values)) || any(decomposition$values < -1e-12)){
    stop("Nominal population logit covariance is invalid.",call.=FALSE)
  }
  active <- which(decomposition$values > .Machine$double.eps *
    max(1,max(decomposition$values)))
  rank <- length(active)
  bytes <- 8*(12*32^rank*max(1,rank) + 12*as.double(nrow(eta))*(d+1))
  if(!is.finite(bytes) ||
     bytes > memoryMB*1048576){
    stop("Nominal population integration exceeds execution$denseMemoryMB.",
      call.=FALSE)
  }
  if(!length(active)){
    result <- .laplace_nominal_probabilities(eta,specification)
  }else{
    root <- decomposition$vectors[,active,drop=FALSE] %*%
      diag(sqrt(decomposition$values[active]),nrow=rank)
    integrate <- function(order){
      rule <- .laplace_normal_quadrature(order)
      grid <- as.matrix(expand.grid(rep(list(seq_len(order)),rank)))
      weights <- apply(grid,1L,function(index) prod(rule$weights[index]))
      shifts <- matrix(rule$nodes[grid],nrow=nrow(grid)) %*% t(root)
      result <- matrix(0,nrow(eta),d+1L)
      for(row in seq_len(nrow(grid))){
        probability <- .laplace_nominal_probabilities(
          sweep(eta,2L,shifts[row,],"+"),specification)
        result <- result+weights[row]*probability
      }
      result
    }
    coarse <- integrate(24L)
    result <- integrate(32L)
    error <- max(abs(result-coarse))
    if(!is.finite(error) || error > 1e-6){
      stop("Nominal population quadrature did not meet the 1e-6 absolute probability accuracy check; use a different integration method.",
        call.=FALSE)
    }
  }
  result
}

.laplace_nominal_control <- function(control){
  if(!inherits(control, "sommer_mmes_control")){
    stop("Nominal preparation requires a control created by mmes_control().",
      call.=FALSE)
  }
  supplied <- attr(control, "supplied")
  native <- attr(control, "native")
  control <- do.call(.mmes_control_build, unclass(control))
  if(!is.null(native)) attr(control, "native") <- native
  if(is.null(supplied)) supplied <- c("formulation", .mme_control_paths(control))
  .mmes_check_compatibility(control, supplied)
  attr(control, "supplied") <- supplied
  .mmes_check_family_controls(control, supplied, stats::gaussian(), NULL)
  if(control$task$glmm$method != "laplace" || control$task$objective != "ml" ||
     control$formulation != "henderson" || control$task$mode != "fit"){
    stop("Internal nominal preparation requires Henderson native Laplace ML in fit mode.",
      call.=FALSE)
  }
  if(control$output$inverse != 0L){
    stop("Nominal uncertainty is post-fit only; set output$inverse=0.",
      call.=FALSE)
  }
  if(length(control$task$glmm$familyParameters)){
    stop("Nominal reference coding has no dispersion or cutpoint parameters; familyParameters must be empty.",
      call.=FALSE)
  }
  control
}

.laplace_prepare_nominal <- function(fixed, random, data=NULL, reference=NULL,
    control=mmes_control(task=list(objective="ml", glmm=list(method="laplace"))),
    naMethodY="exclude", naMethodX="exclude", naMethodRandom="exclude",
    contrasts=NULL){
  control <- .laplace_nominal_control(control)
  scope <- .mmes_blas_scope(control$execution$blasThreads)
  on.exit(scope$restore(), add=TRUE)
  if((!is.null(data) && !is.data.frame(data)) ||
     !inherits(fixed, "formula") || length(fixed) != 3L ||
     !inherits(random, "formula") || length(random) != 2L){
    stop("Nominal preparation requires a response formula and a one-sided random formula; supplied data must be a data frame.",
      call.=FALSE)
  }
  for(method in c(naMethodY, naMethodX, naMethodRandom)){
    if(!is.character(method) || length(method) != 1L || is.na(method) ||
       !method %in% c("exclude", "omit", "fail")){
      stop("Nominal preparation supports missing-data methods 'exclude', 'omit' and 'fail'.",
        call.=FALSE)
    }
  }
  reserved <- c(".nominalCategory", ".nominalResponse")
  if(any(reserved %in% names(data)) ||
     any(reserved %in% all.vars(fixed)) ||
     ".nominalResponse" %in% all.vars(random)){
    stop("Nominal preparation reserves .nominalCategory and .nominalResponse; use .nominalCategory only in the random formula.",
      call.=FALSE)
  }
  frame <- tryCatch(stats::model.frame(fixed, data=data,
    na.action=stats::na.pass, drop.unused.levels=FALSE),
    error=function(condition){
      stop("Nominal response and fixed predictors must be available in the supplied data frame or formula environment: ",
        conditionMessage(condition), call.=FALSE)
    })
  if(!is.null(data) && nrow(frame) != nrow(data)){
    stop("Nominal formula and data must describe the same observations.", call.=FALSE)
  }
  n <- nrow(frame)
  if(is.null(data)){
    data <- as.data.frame(frame)
    environment <- .mmes_formula_env(random, .mmes_formula_env(fixed, parent.frame()))
    for(variable in setdiff(all.vars(random), names(data))){
      value <- tryCatch(get(variable, envir=environment, inherits=TRUE),
        error=function(condition) NULL)
      if(is.null(value)) next
      if((is.atomic(value) || is.factor(value)) && is.null(dim(value)) &&
         length(value) == n){
        data[[variable]] <- value
      }else if(is.matrix(value) && nrow(value) == n){
        data[[variable]] <- I(value)
      }
    }
  }
  nominal <- .laplace_nominal_response(frame[[1L]], reference)
  d <- length(nominal$predictorLevels)
  memoryBytes <- control$execution$denseMemoryMB * 1048576
  if(as.double(n) * d > .Machine$integer.max ||
     2 * d * as.double(utils::object.size(data)) + 8 * n * d * 16 > memoryBytes){
    stop("Nominal stacked-data preparation exceeds execution$denseMemoryMB or sparse index limits.",
      call.=FALSE)
  }
  stackedRows <- rep(seq_len(n), d)
  stacked <- data[stackedRows, , drop=FALSE]
  stacked$.nominalCategory <- factor(rep(nominal$predictorLevels, each=n),
    levels=nominal$predictorLevels)
  stacked$.nominalResponse <- nominal$response[stackedRows]
  # Expand '.' before adding internal columns, so they cannot become covariates.
  stackedFixed <- stats::formula(stats::terms(fixed, data=data))
  stackedFixed[[2L]] <- as.name(".nominalResponse")
  prepared <- .mmes_prepare_model(stackedFixed, random, NULL, stacked,
    stats::gaussian(), naMethodY,
    list(naMethodX=naMethodX, naMethodRandom=naMethodRandom, naMethodR="fail",
      contrasts=contrasts, engine=attr(control, "native")))
  retained <- matrix(prepared$obsInfo$included, n, d)
  if(any(retained != retained[, 1L])){
    stop("Nominal missing-data filtering must retain or exclude every category row of an observation together.",
      call.=FALSE)
  }
  included <- retained[, 1L]
  counts <- tabulate(nominal$response[included], nbins=d + 1L)
  if(any(counts == 0L)){
    stop("Nominal preparation requires retained observations in every declared category.",
      call.=FALSE)
  }
  retainedFrame <- droplevels(frame[included, , drop=FALSE])
  fixedTerms <- stats::delete.response(attr(frame, "terms"))
  basis <- stats::model.matrix(fixedTerms, retainedFrame, contrasts.arg=contrasts)
  if(any(!is.finite(basis))){
    stop("Nominal preparation requires finite retained fixed predictors.", call.=FALSE)
  }
  if(ncol(basis) && Matrix::rankMatrix(basis)[1L] != ncol(basis)){
    stop("Nominal preparation requires a full-rank fixed-effects design.", call.=FALSE)
  }
  if(8 * as.double(nrow(basis)) * max(1, ncol(basis)) * (2 * d + 2) > memoryBytes){
    stop("Nominal fixed design exceeds execution$denseMemoryMB.", call.=FALSE)
  }
  design <- Matrix::kronecker(Matrix::Diagonal(d),
    Matrix::Matrix(basis, sparse=TRUE))
  columns <- data.frame(category=rep(nominal$predictorLevels, each=ncol(basis)),
    coefficient=rep(colnames(basis), d), assignment=rep(attr(basis, "assign"), d),
    stringsAsFactors=FALSE)
  colnames(design) <- make.unique(paste(columns$category, columns$coefficient, sep=":"))
  randomCount <- length(prepared$rtermss)
  if(!randomCount){
    stop("Nominal native Laplace preparation requires random effects.", call.=FALSE)
  }
  randomDesign <- Matrix::Matrix(do.call(cbind, prepared$Z), sparse=TRUE)
  if(any(!is.finite(randomDesign))){
    stop("Nominal preparation requires finite retained random designs.", call.=FALSE)
  }
  offset <- stats::model.offset(frame)
  if(is.null(offset)) offset <- numeric(n)
  offset <- matrix(rep(offset[included], d), sum(included), d,
    dimnames=list(NULL, nominal$predictorLevels))
  if(any(!is.finite(offset))){
    stop("Nominal preparation requires finite retained offsets.", call.=FALSE)
  }
  nominal$response <- nominal$response[included]
  nominal$counts <- stats::setNames(counts, nominal$nativeLevels)
  nominal$fixed <- list(formula=fixed, terms=fixedTerms, columns=colnames(basis),
    contrasts=attr(basis, "contrasts"),
    xlevels=stats::.getXlevels(fixedTerms, retainedFrame), coefficientMap=columns)
  nominal$random <- list(formula=random, terms=prepared$rtermss,
    blockColumns=lapply(prepared$Z, colnames), blockTerms=prepared$Zind)
  nominal$rowMap <- data.frame(originalRow=rep(which(included), d),
    category=rep(nominal$predictorLevels, each=sum(included)),
    stringsAsFactors=FALSE)
  nominal$trainingData <- data[included, , drop=FALSE]
  randomVariables <- setdiff(unique(all.vars(random)), ".nominalCategory")
  randomFactorLevels <- lapply(randomVariables, function(variable){
    value <- data[[variable]]
    if(is.factor(value)) levels(value) else if(is.character(value))
      unique(value[!is.na(value)]) else NULL
  })
  names(randomFactorLevels) <- randomVariables
  randomFactorLevels <- Filter(Negate(is.null), randomFactorLevels)
  randomBlockWidths <- vapply(seq_len(randomCount), function(termIndex)
    sum(vapply(prepared$Z[prepared$Zind == termIndex], ncol, integer(1))),
    integer(1))
  structures <- prepared$covStruct[seq_len(randomCount)]
  parameterNames <- unlist(Map(function(term, structure)
    if(any(structure$free)) paste(term, structure$par_names[structure$free], sep=":")
    else character(),
    prepared$rtermss, structures), use.names=FALSE)
  covarianceStart <- unlist(lapply(structures, function(structure)
    structure$par[structure$free]), use.names=FALSE)
  initial <- c(stats::setNames(numeric(ncol(design)), colnames(design)),
    stats::setNames(covarianceStart, parameterNames))
  structure(list(nominal=nominal, fixedDesign=design, randomDesign=randomDesign,
    offset=offset, covStruct=structures, relationships=prepared$Ai,
    initialParameters=initial,
    obsInfo=prepared$obsInfo[seq_len(n), , drop=FALSE],
    retainedFrame=retainedFrame, randomFactorLevels=randomFactorLevels,
    randomBlockWidths=randomBlockWidths,
    randomDesignBlockWidths=vapply(prepared$Z, ncol, integer(1)),
    control=control),
    class="sommer_nominal_preparation")
}

.laplace_nominal_objective <- function(prepared){
  if(!inherits(prepared, "sommer_nominal_preparation")){
    stop("Nominal objectives require an internal nominal preparation.", call.=FALSE)
  }
  control <- .laplace_nominal_control(prepared$control)
  scope <- .mmes_blas_scope(control$execution$blasThreads)
  on.exit(scope$restore(), add=TRUE)
  context <- sommer_multinomial_objective_create(prepared$nominal$response,
    prepared$fixedDesign, prepared$randomDesign, prepared$offset,
    prepared$covStruct, prepared$relationships, control$execution$denseMemoryMB,
    backend=control$task$glmm$backend,warmStart=control$task$glmm$warmStart,
    trace=control$task$glmm$trace,blocks=control$task$glmm$blocks)
  list(context=context, initialParameters=prepared$initialParameters,
    evaluate=function(parameters, gradient=TRUE){
      gradient <- .mme_flag(gradient, "gradient")
      if(!is.numeric(parameters) || is.complex(parameters) || !is.null(dim(parameters)) ||
         length(parameters) != length(prepared$initialParameters)){
        stop("Nominal objective requires one vector of fixed and free covariance coordinates.",
          call.=FALSE)
      }
      scope <- .mmes_blas_scope(control$execution$blasThreads)
      on.exit(scope$restore(), add=TRUE)
      result <- .laplace_constrained_evaluate(parameters,
        prepared$parameterLayout, function(expanded)
          sommer_multinomial_objective_eval(context, expanded,
            control$task$glmm$modeMaxit, control$task$glmm$modeTol,
            gradient && control$task$glmm$gradient=="analytic"))
      if(gradient && control$task$glmm$gradient=="numeric" && isTRUE(result$valid)){
        score <- .laplace_numeric_gradient(parameters,function(proposal)
          .laplace_constrained_evaluate(proposal,prepared$parameterLayout,
            function(expanded) sommer_multinomial_objective_eval(context,expanded,
              control$task$glmm$modeMaxit,control$task$glmm$modeTol,FALSE))$logLik)
        result <- .laplace_constrained_evaluate(parameters,prepared$parameterLayout,
          function(expanded) sommer_multinomial_objective_eval(context,expanded,
            control$task$glmm$modeMaxit,control$task$glmm$modeTol,FALSE))
        if(!isTRUE(result$valid))
          stop("Numerical nominal score could not restore its accepted mode: ",
            result$reason,".",call.=FALSE)
        result$gradient <- score
      }
      if(!is.null(result$gradient)){
        result$gradient <- stats::setNames(as.numeric(result$gradient),
          names(prepared$initialParameters))
      }
      result
    })
}

.laplace_nominal_diagnostics <- function(optimizer, gradient, boundaries, probability,
    structures=list()){
  scoreMax <- if(length(gradient)) max(abs(gradient)) else 0
  saturated <- sum(apply(probability, 1L, max) >= 1-sqrt(.Machine$double.eps))
  issues <- character()
  if(optimizer$convergence != 0L){
    issues <- c(issues, paste0("Optimizer did not converge (", optimizer$message,
      "); inspect scoreMax and covariance boundaries before using inference."))
  }
  if(length(boundaries)){
    issues <- c(issues, paste0("Covariance boundary: ", paste(boundaries, collapse=", "),
      ". Standard errors and PEV are unavailable; consider a simpler random covariance model."))
  }
  if(saturated){
    issues <- c(issues, paste0(saturated,
      " near-deterministic fitted rows; possible separation or weak identification, not a proof of separation."))
  }
  covarianceCondition <- vapply(structures, function(structure){
    values <- eigen(evaluate_covstruct_cpp(structure, structure$par),
      symmetric=TRUE, only.values=TRUE)$values
    max(values)/min(values)
  }, numeric(1))
  if(any(covarianceCondition > 1/.Machine$double.eps^(1/3))){
    issues <- c(issues, "Ill-conditioned random covariance; inference may be sensitive even if the optimizer converged. Consider a simpler covariance model.")
  }
  list(optimizerCode=optimizer$convergence, optimizerMessage=optimizer$message,
    scoreMax=scoreMax, saturatedRows=saturated, issues=issues,
    covarianceCondition=covarianceCondition,
    uncertaintyStatus=if(optimizer$convergence != 0L) "unavailable: optimizer did not converge"
      else if(length(boundaries)) "unavailable: covariance boundary"
      else "not calculated; use vcov(), summary() or postPEV()")
}

.mmes_laplace_nominal <- function(fixed, random, rcov, data, family, W, weights, vcc,
      REML, maxit, modeMaxit, modeTol, naMethodY, henderson, solveOnly,
      computeCi, mmesArgs, setup=FALSE,backend="auto",warmStart="off",
      gradientMethod="analytic",trace="panels",blocks="off"){
    if(isTRUE(REML)){
      stop("Native multinomial Laplace fits use maximum likelihood; set task$objective='ml'.",
        call.=FALSE)
    }
    if(!isTRUE(henderson) || isTRUE(solveOnly) || as.integer(computeCi) != 0L){
      stop("Native multinomial fitting requires Henderson fit mode and defers uncertainty; set output$inverse=0.",
        call.=FALSE)
    }
    if(!is.null(W) || !is.null(weights)){
      stop("Native multinomial fitting does not support W or weights formulas.",
        call.=FALSE)
    }
    if(!is.null(mmesArgs$naMethodR) &&
       !identical(mmesArgs$naMethodR, "exclude") &&
       !identical(mmesArgs$naMethodR, "omit") &&
       !identical(mmesArgs$naMethodR, "fail")){
      stop("Native multinomial preparation supports naMethodR='exclude', 'omit' or 'fail'.",
        call.=FALSE)
    }
    if(!is.null(mmesArgs$naMethodR) &&
       !identical(mmesArgs$naMethodR, "exclude")){
      stop("Nominal likelihoods have no residual observation coordinates; use naMethodR='exclude'.",
        call.=FALSE)
    }
    if(!is.null(rcov) &&
       !identical(attr(stats::terms(rcov), "term.labels"), "units")){
      stop("Nominal likelihoods use independent observation errors; structured rcov terms are unsupported.",
        call.=FALSE)
    }
    if(is.null(random)){
      stop("Native multinomial Laplace fitting requires at least one random-effect term.",
        call.=FALSE)
    }
    if(identical(mmesArgs$engine$updateMethod, "conditional-fa") ||
       identical(mmesArgs$engine$conditionalEquations, TRUE) ||
       !identical(mmesArgs$engine$factorScoreAugmentation, "none") ||
       !identical(mmesArgs$engine$compression, "off") ||
       isTRUE(mmesArgs$engine$gpu) ||
       !identical(mmesArgs$engine$solver, "auto") ||
       !identical(mmesArgs$engine$acceleration, "none") ||
       !identical(mmesArgs$engine$updateMethod, "ai") ||
       !identical(mmesArgs$engine$fixedRank, "error") ||
       !identical(mmesArgs$engine$covarianceInformation, "full") ||
       !isTRUE(mmesArgs$engine$retainCoefficientMatrix) ||
       any(vapply(mmesArgs$engine$experimental, function(mode)
         identical(mode, "force"), logical(1)))){
      stop("Native multinomial Laplace fitting does not use non-default Gaussian solver, acceleration, update, fixed-rank, conditional FA, factor-score, compression, GPU, or forced experimental controls.",
        call.=FALSE)
    }
    control <- mmes_control(task=list(objective="ml",
        glmm=list(method="laplace", maxit=maxit,
          modeMaxit=modeMaxit, modeTol=modeTol,backend=backend,warmStart=warmStart,
          gradient=gradientMethod,trace=trace,blocks=blocks)),
      execution=list(denseMemoryMB=mmesArgs$engine$denseMemoryMB,
        blasThreads=mmesArgs$blasThreads),
      diagnostics=list(verbose=FALSE, dateWarning=FALSE))
    attr(control, "native") <- mmesArgs$engine
    prepared <- .laplace_prepare_nominal(fixed, random, data, family$reference,
      control=control, naMethodY=naMethodY,
      naMethodX=mmesArgs$naMethodX, naMethodRandom=mmesArgs$naMethodRandom,
      contrasts=mmesArgs$contrasts)
    mapping <- .laplace_covariance_map(prepared$covStruct,
      prepared$nominal$random$terms, vcc)
    prepared$covStruct <- mapping$structures
    freePositions <- lapply(prepared$covStruct, function(structure) which(structure$free))
    parameterStarts <- c(0L, head(cumsum(lengths(freePositions)), -1L))
    fixedCount <- ncol(prepared$fixedDesign)
    covarianceCount <- ncol(mapping$map$T)
    prepared$parameterLayout <- list(fixedCount=fixedCount,
      covarianceCount=covarianceCount, freePositions=freePositions,
      parameterStarts=parameterStarts, freeFamilies=integer(),
      covarianceMap=mapping$map)
    prepared$initialParameters <- c(prepared$initialParameters[seq_len(fixedCount)],
      mapping$map$initial)
    if(setup){
      prepared$vcParams <- mapping$table
      return(prepared)
    }
    objective <- .laplace_nominal_objective(prepared)
    initial <- objective$initialParameters
    initialEvaluation <- objective$evaluate(initial, gradient=FALSE)
    if(!isTRUE(initialEvaluation$valid)){
      stop("Native multinomial starting parameters are invalid: ",
        initialEvaluation$reason, ".", call.=FALSE)
    }
    rejected <- integer()
    value <- function(parameters){
      evaluation <- objective$evaluate(parameters, gradient=FALSE)
      if(!isTRUE(evaluation$valid)){
        reason <- evaluation$reason
        if(!reason %in% names(rejected)) rejected[reason] <<- 0L
        rejected[reason] <<- rejected[reason] + 1L
        return(Inf)
      }
      -evaluation$logLik
    }
    gradient <- function(parameters){
      evaluation <- objective$evaluate(parameters, gradient=TRUE)
      if(!isTRUE(evaluation$valid)){
        stop("Native multinomial score requested at invalid parameters: ",
          evaluation$reason, ".", call.=FALSE)
      }
      -as.numeric(evaluation$gradient)
    }
    optimization <- if(length(initial)){
      stats::nlminb(initial, value, gradient,
        control=list(iter.max=maxit, eval.max=max(1000, 10 * as.double(maxit))))
    }else list(par=numeric(), objective=-initialEvaluation$logLik,
      convergence=0L, message="no estimated parameters")
    final <- objective$evaluate(optimization$par, gradient=TRUE)
    if(!isTRUE(final$valid)){
      stop("Native multinomial optimization ended at invalid parameters: ",
        final$reason, ".", call.=FALSE)
    }
    converged <- identical(optimization$convergence, 0L)
    if(!converged){
      warning("Native multinomial optimizer did not converge: ",
        if(is.null(optimization$message)) paste("code", optimization$convergence)
        else optimization$message, call.=FALSE)
    }
    fixedCount <- ncol(prepared$fixedDesign)
    covarianceCount <- length(initial) - fixedCount
    fittedStructures <- .laplace_proposal_structures(prepared$covStruct,
      optimization$par, prepared$parameterLayout)
    fixedEffects <- stats::setNames(optimization$par[seq_len(fixedCount)],
      colnames(prepared$fixedDesign))
    names(fixedEffects) <- colnames(prepared$fixedDesign)
    nativeProbability <- final$probability
    originalProbability <- nativeProbability[
      , match(prepared$nominal$levels, prepared$nominal$nativeLevels), drop=FALSE]
    colnames(originalProbability) <- prepared$nominal$levels
    rownames(originalProbability) <- rownames(prepared$nominal$trainingData)
    linearPredictors <- matrix(prepared$offset +
      matrix(prepared$fixedDesign %*% fixedEffects,
        nrow=length(prepared$nominal$response), ncol=length(prepared$nominal$predictorLevels)) +
      matrix(prepared$randomDesign %*% final$mode,
          nrow=length(prepared$nominal$response), ncol=length(prepared$nominal$predictorLevels)),
      nrow=length(prepared$nominal$response),
      ncol=length(prepared$nominal$predictorLevels),
      dimnames=list(rownames(prepared$nominal$trainingData),
        prepared$nominal$predictorLevels))
    randomFactorLevels <- prepared$randomFactorLevels
    parameterTable <- .vc_param_table(fittedStructures, prepared$nominal$random$terms)
    covarianceTable <- data.frame(term=parameterTable$term,
      parameter=parameterTable$parameter, estimate=parameterTable$start,
      free=parameterTable$free, StdError=NA_real_, Zratio=NA_real_,
      status="Laplace estimate", stringsAsFactors=FALSE)
    covarianceParameters <- lapply(fittedStructures, function(structure){
      kinds <- lapply(seq_along(structure$par), function(position)
        .vc_kind(structure, position))
      stats::setNames(vapply(seq_along(structure$par), function(position)
        .vc_to_natural(structure$par[position], kinds[[position]]$kind,
          kinds[[position]]$lower, kinds[[position]]$upper), numeric(1)),
        structure$par_names)
    })
    names(covarianceParameters) <- prepared$nominal$random$terms
    boundaries <- .mmes_laplace_boundary(fittedStructures, prepared$covStruct)
    boundaryTable <- matrix(NA_real_, nrow=nrow(covarianceTable),
      ncol=nrow(covarianceTable),
      dimnames=list(paste(covarianceTable$term,covarianceTable$parameter,sep=":"),
        paste(covarianceTable$term,covarianceTable$parameter,sep=":")))
    descriptor <- list(family="multinomial", link="softmax",
      reference=prepared$nominal$reference, levels=prepared$nominal$levels,
      free=FALSE)
    fittedFamily <- multinomial(prepared$nominal$reference)
    result <- list(
      call=match.call(), family=fittedFamily, b=fixedEffects,
      betaSummary=cbind(Estimate=fixedEffects, `Std. Error`=NA_real_,
        `z value`=NA_real_), fitted.values=originalProbability,
      category.probabilities=originalProbability,
      fixedCovariance=NULL, linear.predictors=linearPredictors,
      y=factor(prepared$nominal$originalResponse[prepared$obsInfo$included],
        levels=prepared$nominal$levels),
      offset=prepared$offset,
      AIC=-2*final$logLik+2*length(optimization$par),
      BIC=-2*final$logLik+log(length(prepared$nominal$response))*length(optimization$par),
      llik=as.numeric(final$logLik), deviance=sum(.laplace_nominal_deviance(
        linearPredictors, prepared$nominal, prepared$nominal$response)),
      laplaceFamilies=list(fittedFamily), laplaceFamilyIndex=rep(1L,
        length(prepared$nominal$response)), familySpecifications=list(descriptor),
      familyParameters=NULL, zeroCoefficients=NULL, priorWeights=rep(1,
        length(prepared$nominal$response)), dispersion=1,
      convergence=converged, covParNative=covarianceTable,
      covParCovariance=boundaryTable, engine="native Laplace", REML=FALSE,
      nobs=length(prepared$nominal$response), npar=length(optimization$par),
      fixedFormula=fixed, fixedTerms=prepared$nominal$fixed$terms,
      fixedXlevels=prepared$nominal$fixed$xlevels,
      predictionColumns=prepared$nominal$fixed$columns,
      fixedAssignment=prepared$nominal$fixed$coefficientMap$assignment,
      randomFormula=random,       randomBlockColumns=prepared$nominal$random$blockColumns,
      randomBlockWidths=prepared$randomBlockWidths,
      randomDesignBlockWidths=prepared$randomDesignBlockWidths,
      randomFactorLevels=randomFactorLevels,
      predictionSettings=list(nominal=TRUE, mmesArgs=mmesArgs,
        naMethodY=naMethodY),
      fixedDesign=prepared$fixedDesign, randomDesign=prepared$randomDesign,
      randomMode=as.numeric(final$mode), laplaceNativeResponse=prepared$nominal$response,
      laplaceModeHessian=final$modeHessian,
      laplaceObservationCurvature=final$observationCurvature,
      randomMetadata=prepared$nominal$random$terms,
      covStruct=fittedStructures, randomRelationships=prepared$relationships,
      covPar=covarianceParameters, u=matrix(as.numeric(final$mode),ncol=1L),
      contrasts=prepared$nominal$fixed$contrasts, optimizer=optimization,
      laplaceParameterLayout=prepared$parameterLayout,
      vcParams=mapping$table, vcc=vcc,
      obsInfo=prepared$obsInfo, nominal=prepared$nominal,
      nominalPreparation=prepared,
      laplaceDiagnostics=list(outerConverged=converged, modeConverged=TRUE,
        modeScore=final$modeScore, curvatureValid=NA,
        curvatureReason="not calculated; use vcov() or summary()",
        uncertaintyValid=FALSE, uncertaintyComputed=FALSE,
        boundary=length(boundaries)>0L, boundaryParameters=boundaries,
        initialGLMConverged=NA, gradient=as.numeric(final$gradient),
        nativeWork=final$work, extremeProbabilities=any(
          apply(originalProbability,1,max)>1-.Machine$double.eps),
        rejectedEvaluations=rejected),
      modeIterations=final$iterations, logDetHessian=final$logDetHessian,
      executionPlan=list(glmm=list(selected="Native multinomial Laplace")),
      laplacePostPEVSettings=list(inverse=0L, panelSize=as.integer(
        mmesArgs$engine$experimental$panelCi2),
        denseMemoryMB=mmesArgs$engine$denseMemoryMB, modeMaxit=modeMaxit,
        modeTol=modeTol,warmStart=warmStart,gradient=gradientMethod,trace=trace,blocks=blocks))
    result$laplaceDiagnostics$boundaryParameters <- unique(boundaries)
    result$laplaceDiagnostics <- c(result$laplaceDiagnostics,
      .laplace_nominal_diagnostics(optimization, final$gradient, unique(boundaries),
        originalProbability, fittedStructures))
    result$nominal$specification <- descriptor
    result$nominal$counts <- stats::setNames(prepared$nominal$counts,
      prepared$nominal$nativeLevels)
    result$randomBlockColumns <- lapply(seq_along(prepared$nominal$random$blockColumns),
      function(position) prepared$nominal$random$blockColumns[[position]])
    class(result) <- c("mmes.glmm.laplace", "mmes.glmm", "mmes")
    result$Ci <- result$VarU <- result$uPevList <- NULL
    result$uPevMode <- 0L
    result$nominalPreparation$control$task$glmm$backend <- final$work$backend$selected
    if(result$laplaceDiagnostics$extremeProbabilities){
      warning("Native multinomial fitted probabilities are numerically 0 or 1; check separation and identifiability.",
        call.=FALSE)
    }
    result
}

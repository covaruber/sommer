.nonlinear_prepare <- function(fixed, random, data, control, naMethodX,
    naMethodY, naMethodRandom, naMethodR, contrasts, input=NULL){
  args <- c(list(fixed=fixed, random=random, data=data,
    naMethodX=naMethodX, naMethodY=naMethodY,
    naMethodRandom=naMethodRandom, naMethodR=naMethodR, contrasts=contrasts),
    .mmes_control_args(control),
    list(.prepareOnly=TRUE, .internal=list(nonlinearInput=input)))
  args$returnParam <- TRUE
  args$verbose <- FALSE
  do.call(.mmes_fit, args)
}

.nonlinear_identification <- function(fixed, jacobian){
  design <- cbind(fixed, Matrix::Matrix(jacobian, sparse=TRUE))
  basis <- fixed_design_basis_cpp(
    as(as(as(design, "dMatrix"), "generalMatrix"), "CsparseMatrix"))
  !length(basis$omitted)
}

.nonlinear_covariance_parameters <- function(structures){
  lapply(structures, function(structure){
    vapply(seq_along(structure$par), function(position){
      kind <- .vc_kind(structure, position)
      .vc_to_natural(structure$par[position], kind$kind, kind$lower, kind$upper)
    }, numeric(1))
  })
}

.nonlinear_graph_start <- function(graph,settings){
  start <- settings$start
  if(is.null(start))
    stop("Nonlinear graph fitting requires a complete named task$nonlinear$start vector.",
      call.=FALSE)
  if(!is.numeric(start) || is.complex(start) || !is.null(dim(start)) ||
     length(start)!=length(graph$parameters) || is.null(names(start)) ||
     anyNA(names(start)) || anyDuplicated(names(start)) ||
     !setequal(names(start),graph$parameters) || any(!is.finite(start))){
    stop("task$nonlinear$start must be a finite named vector containing exactly: ",
      paste(graph$parameters,collapse=", "),".",call.=FALSE)
  }
  stats::setNames(as.numeric(start[graph$parameters]),graph$parameters)
}

.mmes_nonlinear_fit <- function(fixed, random=NULL, data, nonlinear, control,
    rcov=NULL, W=NULL, weights=NULL, family=stats::gaussian(), vcc=NULL,
    covPar=NULL, naMethodX="exclude", naMethodY="exclude",
    naMethodRandom="exclude", naMethodR="exclude", contrasts=NULL){
  if(!.nonlinear_exact_gaussian(family, control)){
    return(.mmes_nonlinear_family_fit(fixed, random, data, nonlinear, control,
      family, vcc, naMethodX, naMethodY, naMethodRandom, naMethodR, contrasts))
  }
  profileStart <- if(identical(control$diagnostics$profile,"force")) proc.time()[["elapsed"]] else NULL
  graphFit <- inherits(nonlinear,"sommer_nonlinear_graph")
  model <- if(graphFit) .nonlinear_graph_validate(nonlinear) else
    .nonlinear_validate(nonlinear)
  if(!graphFit) random <- .nonlinear_joint_formula(random, model, data)
  graphLayout <- if(graphFit) .nonlinear_graph_native_layout(model,data,
    allowMissing=TRUE) else NULL
  input <- if(graphFit) graphLayout$inputs else
    .nonlinear_input(model, data, missingOK=TRUE)
  scope <- .mmes_blas_scope(control$execution$blasThreads)
  on.exit(scope$restore(), add=TRUE)
  prepared <- .nonlinear_prepare(fixed, random, data, control, naMethodX,
    naMethodY, naMethodRandom, naMethodR, contrasts, input)
  included <- prepared$obsInfo$included
  x <- prepared$X
  z <- Matrix::Matrix(do.call(cbind, prepared$Z), sparse=TRUE)
  curveRandom <- if(graphFit) NULL else .nonlinear_random_layout(model, prepared)
  prepared$rtermss <- .nonlinear_term_labels(model,prepared$rtermss)
  y <- stats::model.response(prepared$modelFrame)
  if(!is.numeric(y) || !is.null(dim(y))){
    stop("Nonlinear fitting requires one numeric response.", call.=FALSE)
  }
  y <- as.numeric(y[included])
  offset <- stats::model.offset(prepared$modelFrame)
  offset <- if(is.null(offset)) numeric(length(y)) else offset[included]
  if(!graphFit) input <- input[included]
  retainedInput <- if(graphFit) input[included,,drop=FALSE] else input
  if(any(!is.finite(y)) || any(!is.finite(offset)) ||
     any(!is.finite(retainedInput)) ||
     any(!is.finite(x)) || any(!is.finite(z))){
    stop("Nonlinear fitting requires finite retained responses, inputs, offsets and designs.",
      call.=FALSE)
  }
  k <- length(model$parameters)
  p <- ncol(x)
  if(length(y) <= p+k || (!graphFit && length(unique(input)) < k)){
    stop("Too few distinct inputs or observations to identify the nonlinear mean and residual variance.",
      call.=FALSE)
  }
  if(8*(40*as.double(length(y))+12*as.double(ncol(z))+
        as.double(length(y))*(p+k)) > control$execution$denseMemoryMB*1048576){
    stop("Nonlinear preparation workspace exceeds execution$denseMemoryMB.", call.=FALSE)
  }
  theta <- if(graphFit) .nonlinear_graph_start(model,control$task$nonlinear) else
    .nonlinear_start(model, input, y-offset, control$task$nonlinear)
  initialPredictor <- if(graphFit){
    .nonlinear_graph_first(model,data[included,,drop=FALSE],theta,control)
  }else sommer_nonlinear_curve_eval(model$curve,input,theta)
  initialMean <- if(graphFit) initialPredictor$output else initialPredictor$mean
  initialJacobian <- initialPredictor$jacobian
  if(!.nonlinear_identification(x, initialJacobian)){
    stop("The combined fixed/nonlinear mean is locally rank deficient at the start; revise the fixed formula or supply identifiable nonlinear starts.",
      call.=FALSE)
  }
  beta <- if(p) stats::lm.fit(as.matrix(x),
    y-offset-as.numeric(initialMean))$coefficients else numeric()
  initialResidual <- y-offset-as.numeric(initialMean)-as.numeric(x %*% beta)
  variance <- mean(initialResidual^2)
  if(!is.finite(variance) || variance <= 0){
    stop("Nonlinear initialization produced non-positive residual variance.", call.=FALSE)
  }
  count <- length(prepared$rtermss)
  mapping <- .laplace_covariance_map(prepared$covStruct[seq_len(count)],
    prepared$rtermss, vcc)
  .nonlinear_covariance_eligibility(mapping$structures, prepared$rtermss)
  # Nonlinear predictor coordinates are part of the fixed-mean block, not covariance.
  layout <- list(fixedCount=p+k, covarianceCount=length(mapping$map$initial),
    freePositions=lapply(mapping$structures, function(s) which(s$free)),
    parameterStarts=cumsum(c(0L, head(vapply(mapping$structures,
      function(s) sum(s$free), integer(1)), -1L))),
    covarianceMap=mapping$map)
  parameters <- c(if(p) stats::setNames(beta, paste0("fixed:", colnames(x))),
    stats::setNames(theta,paste0(if(graphFit) "graph:" else "curve:",
      model$parameters)),
    mapping$map$initial, `residual:logVariance`=log(variance))
  if(graphFit){
    context <- sommer_nonlinear_graph_gaussian_create(x,z,y,offset,
      graphLayout$inputs[included,,drop=FALSE],graphLayout$nodes,
      graphLayout$output,graphLayout$transformed,mapping$structures,
      prepared$Ai,control$execution$denseMemoryMB)
    nativeEvaluate <- sommer_nonlinear_gaussian_eval
  }else if(is.null(curveRandom)){
    context <- sommer_nonlinear_gaussian_create(x, z, y, offset, input,
      model$curve, mapping$structures, prepared$Ai, control$execution$denseMemoryMB)
    nativeEvaluate <- sommer_nonlinear_gaussian_eval
  }else if(identical(curveRandom$version,2L)){
    mode <- control$task$nonlinear$mode
    specifications <- .laplace_family_specs(stats::gaussian(),list())$specifications
    context <- sommer_nonlinear_family_create(x,z,y,offset,input,model$curve,
      curveRandom$coordinate[1L],curveRandom$start,mapping$structures,prepared$Ai,
      specifications,rep(0L,length(y)),rep(1,length(y)),
      control$execution$denseMemoryMB,mode$maxIters,mode$scoreTol,mode$stepTol,
      coordinateMap=curveRandom$coordinateMap*1.0)
    nativeEvaluate <- sommer_nonlinear_family_eval
  }else{
    mode <- control$task$nonlinear$mode
    context <- sommer_nonlinear_laplace_create(x, z, y, offset, input, model$curve,
      curveRandom$coordinate, curveRandom$start, mapping$structures, prepared$Ai,
      control$execution$denseMemoryMB, mode$maxIters, mode$scoreTol, mode$stepTol)
    nativeEvaluate <- sommer_nonlinear_laplace_eval
  }
  evaluate <- function(parameters, gradient=TRUE){
    .laplace_constrained_evaluate(parameters, layout, function(expanded)
      nativeEvaluate(context, expanded, gradient))
  }
  preparationSeconds <- if(!is.null(profileStart)) proc.time()[["elapsed"]]-profileStart else NULL
  optimized <- .nonlinear_optimize(parameters,evaluate,control,"Nonlinear Gaussian ML")
  optimization <- optimized$optimization
  final <- optimized$final
  rejected <- optimized$rejected
  converged <- identical(optimization$convergence, 0L)
  theta <- optimization$par[p+seq_len(k)]
  names(theta) <- model$parameters
  curve <- if(graphFit) NULL else sommer_nonlinear_curve_eval(model$curve,input,theta)
  finalJacobian <- if(graphFit)
    .nonlinear_graph_first(model,data[included,,drop=FALSE],theta,control)$jacobian else
      curve$jacobian
  identified <- .nonlinear_identification(x,finalJacobian)
  if(!identified){
    warning("The fitted nonlinear mean is locally rank deficient; curve parameters may not be identifiable.",
      call.=FALSE)
  }
  structures <- .laplace_proposal_structures(mapping$structures,
    optimization$par, layout)
  names(structures) <- prepared$rtermss
  table <- .vc_param_table(structures, prepared$rtermss)
  covariance <- .nonlinear_covariance_parameters(structures)
  names(covariance) <- prepared$rtermss
  beta <- stats::setNames(optimization$par[seq_len(p)], colnames(x))
  if(!is.null(curveRandom)){
    final$fixedMean <- offset+as.numeric(x %*% beta)+as.numeric(curve$mean)
  }
  fixedTerms <- stats::delete.response(attr(prepared$modelFrame, "terms"))
  retainedFrame <- droplevels(prepared$modelFrame[included, , drop=FALSE])
  training <- data[included, , drop=FALSE]
  randomVariables <- intersect(all.vars(random), names(training))
  randomLevels <- lapply(training[randomVariables], function(value){
    if(is.factor(value)) levels(value) else if(is.character(value)) unique(value)
    else NULL
  })
  randomLevels <- Filter(Negate(is.null), randomLevels)
  boundary <- .mmes_laplace_boundary(structures, mapping$structures)
  nonlinearParameters <- if(graphFit){
    natural <- theta
    logCoordinates <- which(vapply(model$parameterMap,function(map)
      identical(map$transform,"log"),logical(1)))
    if(length(logCoordinates)) natural[logCoordinates] <- exp(natural[logCoordinates])
    natural
  }else .nonlinear_natural(model,theta)
  result <- list(
    nonlinear=model, nonlinearGraph=graphFit,
    nonlinearParameters=nonlinearParameters,
    nonlinearCoordinates=stats::setNames(theta, model$parameters),
    nonlinearRandom=curveRandom,
    groupParameters=.nonlinear_random_parameters(model, theta, curveRandom, final$mode),
    likelihoodMethod=if(is.null(curveRandom)) "exact Gaussian ML" else "Laplace Gaussian ML",
    nonlinearParameterLayout=list(version=1L, fixed=names(beta),
      curve=if(graphFit) character() else model$parameters,
      graph=if(graphFit) model$parameters else character(),
      positive=if(graphFit) which(vapply(model$parameterMap,function(map)
        identical(map$transform,"log"),logical(1))) else model$positive,
      optimizer=names(parameters), covariance=layout),
    b=beta, u=as.numeric(final$mode), randomMode=as.numeric(final$mode),
    latentRandomMode=as.numeric(final$latentMode),
    randomCoordinates=list(
      representation=if(.nonlinear_has_reduced_rank(structures)) "latent reduced rank" else "full",
      fullDimension=length(final$mode), latentDimension=length(final$latentMode)),
    fitted.values=as.numeric(final$fitted), fixedMean=as.numeric(final$fixedMean),
    y=y, offset=offset, llik=as.numeric(final$logLik), REML=FALSE,
    nobs=length(y), npar=length(parameters), dispersion=exp(tail(optimization$par, 1L)),
    family=family, convergence=converged,
    engine=if(is.null(curveRandom)) "native nonlinear Gaussian ML" else
      "native nonlinear Laplace Gaussian ML",
    optimizer=optimization, covStruct=structures, covPar=covariance,
    covParNative=table, randomMetadata=prepared$rtermss,
    randomRelationships=prepared$Ai, fixedDesign=x, randomDesign=z,
    fixedFormula=fixed, randomFormula=.nonlinear_prediction_formula(random, training),
    fixedTerms=fixedTerms,
    fixedVariables=intersect(all.vars(fixedTerms), names(training)),
    fixedXlevels=stats::.getXlevels(fixedTerms, retainedFrame),
    predictionColumns=colnames(x), contrasts=attr(x, "contrasts"),
    randomBlockColumns=lapply(prepared$Z, colnames),
    randomBlockWidths=vapply(seq_along(prepared$rtermss),function(term)
      sum(vapply(prepared$Z[prepared$Zind==term],ncol,integer(1))),integer(1)),
    randomDesignBlockWidths=vapply(prepared$Z,ncol,integer(1)),
    randomFactorLevels=randomLevels, predictionData=training[randomVariables],
    input=if(graphFit) input[included,,drop=FALSE] else input,
    graphNative=if(graphFit) graphLayout[c("nodes","output","transformed")] else NULL,
    graphData=if(graphFit) graphLayout$inputs[included,,drop=FALSE] else NULL,
    trainingData=if(graphFit) training else NULL,
    obsInfo=prepared$obsInfo, vcc=vcc,
    predictionSettings=list(naMethodX=naMethodX, naMethodY=naMethodY,
      naMethodRandom=naMethodRandom, naMethodR=naMethodR, contrasts=contrasts),
    nonlinearDiagnostics=list(gradient=stats::setNames(as.numeric(final$gradient),
      names(parameters)), initialParameters=parameters,
      modeScore=final$modeScore, identified=identified,
      modeIterations=final$modeIterations, failedStarts=final$failedStarts,
      failedStartReasons=final$failedStartReasons,
      competingModes=final$competingModes, logDetHessian=final$logDetHessian,
      boundaryParameters=boundary, rejectedEvaluations=rejected,
      work=final$work, uncertaintyComputed=FALSE),
    uPevMode=0L, retainCoefficientMatrix=FALSE, covarianceInformation="none")
  class(result) <- c("mmes.nonlinear", "mmes")
  if(!is.null(curveRandom)){
    result$curveRandomMode <- stats::setNames(
      result$randomMode[curveRandom$start+seq_len(curveRandom$count)], curveRandom$columns)
    result$additiveRandomMode <- result$randomMode[seq_len(curveRandom$start)]
  }
  .nonlinear_fit_diagnostics(result,optimized,control,profileStart,preparationSeconds)
}

.nonlinear_family_prepare <- function(fixed, random, data, family, input,
    control, naMethodX, naMethodY, naMethodRandom, naMethodR, contrasts){
  settings <- control$task$glmm$familyParameters
  singleCategorical <- !inherits(family, "sommer_familym") &&
    family$family %in% c("ordinal", "multinomial")
  if(singleCategorical){
    reserved <- c(".nonlinearResponse", ".nonlinearTrait")
    if(any(reserved %in% names(data)))
      stop("Nonlinear categorical preparation reserves .nonlinearResponse and .nonlinearTrait.",
        call.=FALSE)
    response <- stats::model.response(stats::model.frame(fixed,data,
      na.action=stats::na.pass))
    if(family$family=="multinomial"){
      .laplace_nominal_response(response,family$reference)
    }else if(!is.ordered(response) && (!is.numeric(response) ||
        is.object(response) || is.complex(response) || !is.null(dim(response)))){
      stop("Ordinal responses must be ordered factors or plain numeric integer category codes.",
        call.=FALSE)
    }
    if(is.factor(response)){
      family$.nonlinearLevels <- levels(response)
      response <- as.numeric(response)
    }
    data$.nonlinearResponse <- response
    data$.nonlinearTrait <- "response"
    fixed[[2L]] <- as.name(".nonlinearResponse")
    family <- familym(response=family, by=".nonlinearTrait")
    settings <- list(response=settings)
    if(".nominalCategory" %in% all.vars(random)){
      random <- .nonlinear_replace_symbol(random, ".nominalCategory", ".familyCoordinate")
    }
  }
  complete <- !is.na(data[[input]]) & .laplace_zero_complete(family,data,nrow(data))
  if(naMethodX=="fail" && any(!complete))
    stop("Missing nonlinear inputs or zero-component predictors with naMethodX='fail'.",
      call.=FALSE)
  args <- list(naMethodX=naMethodX,naMethodRandom=naMethodRandom,
    naMethodR=naMethodR,contrasts=contrasts,engine=attr(control,"native"),
    observationMask=complete)
  prepared <- .mmes_prepare_model(fixed,random,NULL,data,family,naMethodY,args,settings)
  categorical <- prepared$categorical
  if(!is.null(categorical)){
    specifications <- categorical$components
    families <- family$families
    workingData <- categorical$expandedData
  }else{
    model <- .laplace_family_specs(family,settings)
    specifications <- model$specifications
    families <- model$families
    workingData <- data
  }
  included <- prepared$obsInfo$included
  index <- .laplace_family_index(family,workingData,nrow(workingData))[included]
  specifications <- .laplace_zero_designs(families,specifications,
    workingData,included,index)
  response <- stats::model.response(prepared$modelFrame)
  declaredTrials <- if(inherits(family,"sommer_familym") && !is.null(family$trials))
    workingData[[family$trials]][included] else NULL
  response <- .laplace_response(response,included,specifications,index,declaredTrials)
  offset <- stats::model.offset(prepared$modelFrame)
  if(is.null(offset)) offset <- numeric(length(included))
  list(prepared=prepared,family=family,families=families,specifications=specifications,
    index=index,response=response,offset=offset[included],data=workingData,
    input=workingData[[input]][included],random=random,
    singleCategorical=singleCategorical,args=args)
}

.nonlinear_replace_symbol <- function(expression, from, to){
  if(is.symbol(expression)) return(if(as.character(expression)==from) as.name(to) else expression)
  if(is.call(expression)){
    result <- as.call(lapply(as.list(expression),.nonlinear_replace_symbol,from,to))
    if(inherits(expression,"formula")){
      class(result) <- "formula"
      environment(result) <- environment(expression)
    }
    return(result)
  }
  expression
}

.nonlinear_optimize <- function(parameters, evaluate, control, label="Nonlinear ML"){
  profiling <- identical(control$diagnostics$profile,"force")
  started <- if(profiling) proc.time()[["elapsed"]] else NULL
  evaluationSeconds <- 0
  measured <- function(parameters){
    if(!profiling) return(evaluate(parameters))
    start <- proc.time()[["elapsed"]]
    value <- evaluate(parameters)
    evaluationSeconds <<- evaluationSeconds+proc.time()[["elapsed"]]-start
    value
  }
  cachedParameters <- cachedValue <- NULL
  rejected <- integer()
  cached <- function(parameters){
    if(!identical(parameters,cachedParameters)){
      cachedValue <<- measured(parameters)
      cachedParameters <<- parameters
    }
    cachedValue
  }
  objective <- function(parameters){
    value <- cached(parameters)
    if(!isTRUE(value$valid)){
      reason <- value$reason
      if(!reason %in% names(rejected)) rejected[reason] <<- 0L
      rejected[reason] <<- rejected[reason]+1L
      return(Inf)
    }
    -value$logLik
  }
  gradient <- function(parameters){
    value <- cached(parameters)
    if(!isTRUE(value$valid))
      stop("Nonlinear gradient requested at invalid parameters: ",value$reason,".",call.=FALSE)
    -as.numeric(value$gradient)
  }
  initial <- cached(parameters)
  if(!isTRUE(initial$valid))
    stop("Invalid nonlinear starting parameters: ",initial$reason,".",call.=FALSE)
  convergence <- control$convergence
  optimization <- stats::nlminb(parameters,objective,gradient=gradient,
    control=list(iter.max=convergence$maxIters,
      eval.max=max(1000,10*as.double(convergence$maxIters)),
      rel.tol=convergence$likelihoodTol,x.tol=convergence$parameterTol))
  final <- measured(optimization$par)
  if(!isTRUE(final$valid))
    stop("Nonlinear optimization ended at invalid parameters: ",final$reason,".",call.=FALSE)
  if(optimization$convergence!=0L)
    warning(label," optimizer did not converge: ",optimization$message,call.=FALSE)
  if(isTRUE(final$competingModes>0L))
    warning("Distinct conditional modes were found; Laplace uses the best deterministic candidate, not a guaranteed global mode.",
      call.=FALSE)
  list(optimization=optimization,final=final,rejected=rejected,
    profile=if(profiling) c(outerTotal=proc.time()[["elapsed"]]-started,
      marginalEvaluations=evaluationSeconds) else NULL)
}

.mmes_nonlinear_family_fit <- function(fixed, random, data, nonlinear, control,
    family, vcc, naMethodX, naMethodY, naMethodRandom, naMethodR, contrasts){
  profileStart <- if(identical(control$diagnostics$profile,"force")) proc.time()[["elapsed"]] else NULL
  model <- .nonlinear_validate(nonlinear)
  random <- .nonlinear_joint_formula(random,model,data)
  .nonlinear_input(model,data,missingOK=TRUE)
  scope <- .mmes_blas_scope(control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  setup <- .nonlinear_family_prepare(fixed,random,data,family,model$input,control,
    naMethodX,naMethodY,naMethodRandom,naMethodR,contrasts)
  prepared <- setup$prepared
  included <- prepared$obsInfo$included
  training <- setup$data[included,,drop=FALSE]
  x <- prepared$X
  z <- Matrix::Matrix(do.call(cbind,prepared$Z),sparse=TRUE)
  curveRandom <- .nonlinear_random_layout(model,prepared)
  prepared$rtermss <- .nonlinear_term_labels(model,prepared$rtermss)
  input <- setup$input
  p <- ncol(x); k <- length(model$parameters)
  if(any(!is.finite(input)) || any(!is.finite(setup$offset)) ||
     any(!is.finite(x)) || any(!is.finite(z)))
    stop("Nonlinear family preparation requires finite retained inputs, offsets and designs.",
      call.=FALSE)
  if(length(unique(input))<k || length(input)<=p+k)
    stop("Too few distinct inputs or observations to identify the nonlinear mean.",call.=FALSE)
  if(8*as.double(length(input))*(150+p+k)>control$execution$denseMemoryMB*1048576)
    stop("Nonlinear family preparation exceeds execution$denseMemoryMB.",call.=FALSE)
  theta <- .nonlinear_start(model,input,rep(1,length(input)),control$task$nonlinear)
  curve <- sommer_nonlinear_curve_eval(model$curve,input,theta)
  locations <- lapply(seq_along(setup$specifications),function(i){
    spec <- setup$specifications[[i]]
    if(spec$family=="ordinal" && spec$free) as.numeric(setup$index==i) else NULL
  })
  location <- do.call(cbind,Filter(Negate(is.null),locations))
  identified <- .nonlinear_identification(cbind(x,location),curve$jacobian)
  if(!identified)
    stop("The combined fixed/nonlinear mean and free cutpoint locations are locally rank deficient; revise the model or starts.",
      call.=FALSE)
  mapping <- .laplace_covariance_map(prepared$covStruct[seq_along(prepared$rtermss)],
    prepared$rtermss,vcc)
  .nonlinear_covariance_eligibility(mapping$structures,prepared$rtermss)
  layout <- list(fixedCount=p+k,covarianceCount=length(mapping$map$initial),
    covarianceMap=mapping$map,
    freePositions=lapply(mapping$structures,function(s) which(s$free)),
    parameterStarts=cumsum(c(0L,head(vapply(mapping$structures,
      function(s) sum(s$free),integer(1)),-1L))))
  familyStarts <- unlist(lapply(setup$specifications,function(spec){
    if(!spec$free) return(numeric())
    if(spec$family=="ordinal") .laplace_ordinal_coordinates(spec$thresholds)
    else if(is.null(spec$zeroDesign)) spec$logParameter else spec$zeroCoefficients
  }),use.names=FALSE)
  familyNames <- unlist(lapply(seq_along(setup$specifications),function(i){
    spec <- setup$specifications[[i]]
    count <- .laplace_family_count(spec)
    if(!count) return(character())
    paste0("family:",names(setup$specifications)[i],":",
      if(spec$family=="ordinal") paste0("cutpoint",seq_len(count)) else
        if(is.null(spec$zeroDesign)) spec$parameter else colnames(spec$zeroDesign))
  }),use.names=FALSE)
  parameters <- c(if(p) stats::setNames(numeric(p),paste0("fixed:",colnames(x))),
    stats::setNames(theta,paste0("curve:",model$parameters)),mapping$map$initial,
    stats::setNames(familyStarts,familyNames))
  mode <- control$task$nonlinear$mode
  context <- sommer_nonlinear_family_create(x,z,setup$response$nativeResponse,
    setup$offset,input,model$curve,if(is.null(curveRandom)) 0L else curveRandom$coordinate[1L],
    if(is.null(curveRandom)) ncol(z) else curveRandom$start,
    mapping$structures,prepared$Ai,setup$specifications,setup$index-1L,
    setup$response$trials,control$execution$denseMemoryMB,
    mode$maxIters,mode$scoreTol,mode$stepTol,
    coordinateMap=if(identical(curveRandom$version,2L)) curveRandom$coordinateMap*1.0 else NULL)
  evaluate <- function(parameters){
    .laplace_constrained_evaluate(parameters,layout,function(expanded)
      sommer_nonlinear_family_eval(context,expanded))
  }
  preparationSeconds <- if(!is.null(profileStart)) proc.time()[["elapsed"]]-profileStart else NULL
  optimized <- .nonlinear_optimize(parameters,evaluate,control)
  optimization <- optimized$optimization; final <- optimized$final
  theta <- optimization$par[p+seq_len(k)]
  curve <- sommer_nonlinear_curve_eval(model$curve,input,theta)
  identified <- .nonlinear_identification(cbind(x,location),curve$jacobian)
  if(!identified)
    warning("The fitted nonlinear mean is locally rank deficient; curve parameters may not be identifiable.",
      call.=FALSE)
  structures <- .laplace_proposal_structures(mapping$structures,optimization$par,layout)
  names(structures) <- prepared$rtermss
  specifications <- .laplace_proposal_specs(setup$specifications,optimization$par,layout)
  families <- .laplace_update_families(setup$families,specifications,setup$index)
  names(families) <- names(specifications)
  fittedFamily <- if(inherits(setup$family,"sommer_familym")){
    value <- setup$family
    value$families <- families
    value
  }else families[[1L]]
  beta <- stats::setNames(optimization$par[seq_len(p)],colnames(x))
  fixedMean <- setup$offset+as.numeric(x%*%beta)+as.numeric(
    sommer_nonlinear_curve_eval(model$curve,input,theta,gradient=FALSE)$mean)
  fixedTerms <- stats::delete.response(attr(prepared$modelFrame,"terms"))
  retainedFrame <- droplevels(prepared$modelFrame[included,,drop=FALSE])
  randomVariables <- intersect(all.vars(setup$random),names(training))
  randomLevels <- lapply(training[randomVariables],function(value)
    if(is.factor(value)) levels(value) else if(is.character(value)) unique(value) else NULL)
  randomLevels <- Filter(Negate(is.null),randomLevels)
  covariance <- .nonlinear_covariance_parameters(structures)
  names(covariance) <- prepared$rtermss
  dispersion <- vapply(specifications,function(s)
    if(s$family %in% c("Gamma","gaussian")) exp(s$logParameter) else 1,numeric(1))
  result <- list(nonlinear=model,nonlinearParameters=.nonlinear_natural(model,theta),
    nonlinearCoordinates=stats::setNames(theta,model$parameters),nonlinearRandom=curveRandom,
    groupParameters=.nonlinear_random_parameters(model,theta,curveRandom,final$mode),
    nonlinearFamily=TRUE,singleCategorical=setup$singleCategorical,
    likelihoodMethod="Laplace family ML",b=beta,u=as.numeric(final$mode),
    randomMode=as.numeric(final$mode),latentRandomMode=as.numeric(final$latentMode),
    randomCoordinates=list(representation=if(.nonlinear_has_reduced_rank(structures))
      "latent reduced rank" else "full",fullDimension=length(final$mode),
      latentDimension=length(final$latentMode)),
    fitted.values=as.numeric(final$fitted),linear.predictors=as.numeric(final$predictor),
    fixedMean=fixedMean,y=setup$response$response,priorWeights=setup$response$trials,
    nativeResponse=setup$response$nativeResponse,
    offset=setup$offset,llik=as.numeric(final$logLik),REML=FALSE,
    nobs=if(is.null(prepared$categorical)) length(input) else
      length(unique(prepared$categorical$rowMap)),npar=length(parameters),
    dispersion=dispersion,family=if(setup$singleCategorical) families[[1L]] else fittedFamily,
    predictionFamily=fittedFamily,laplaceFamilies=families,
    laplaceFamilyIndex=setup$index,familySpecifications=specifications,
    categorical=prepared$categorical,
    convergence=optimization$convergence==0L,engine="native nonlinear Laplace family ML",
    optimizer=optimization,covStruct=structures,covPar=covariance,
    covParNative=.vc_param_table(structures,prepared$rtermss),
    randomMetadata=prepared$rtermss,randomRelationships=prepared$Ai,
    fixedDesign=x,randomDesign=z,fixedFormula=fixed,
    randomFormula=.nonlinear_prediction_formula(setup$random,training),fixedTerms=fixedTerms,
    fixedVariables=intersect(all.vars(fixedTerms),names(training)),
    fixedXlevels=stats::.getXlevels(fixedTerms,retainedFrame),
    predictionColumns=colnames(x),contrasts=attr(x,"contrasts"),
    randomBlockColumns=lapply(prepared$Z,colnames),
    randomBlockWidths=vapply(seq_along(prepared$rtermss),function(term)
      sum(vapply(prepared$Z[prepared$Zind==term],ncol,integer(1))),integer(1)),
    randomDesignBlockWidths=vapply(prepared$Z,ncol,integer(1)),
    randomFactorLevels=randomLevels,predictionData=training[randomVariables],
    input=input,obsInfo=prepared$obsInfo,vcc=vcc,
    predictionSettings=list(naMethodX=naMethodX,naMethodY=naMethodY,
      naMethodRandom=naMethodRandom,naMethodR=naMethodR,contrasts=contrasts,
      mmesArgs=setup$args),
    laplacePostPEVSettings=list(denseMemoryMB=control$execution$denseMemoryMB),
    nonlinearParameterLayout=list(version=2L,fixed=names(beta),curve=model$parameters,
      positive=model$positive,optimizer=names(parameters),covariance=layout),
    nonlinearDiagnostics=list(gradient=stats::setNames(as.numeric(final$gradient),names(parameters)),
      initialParameters=parameters,modeScore=final$modeScore,identified=identified,
      modeIterations=final$modeIterations,failedStarts=final$failedStarts,
      failedStartReasons=final$failedStartReasons,competingModes=final$competingModes,
      logDetHessian=final$logDetHessian,rejectedEvaluations=optimized$rejected,
      boundaryParameters=.mmes_laplace_boundary(structures,mapping$structures),
      work=final$work,uncertaintyComputed=FALSE),
    uPevMode=0L,retainCoefficientMatrix=FALSE,covarianceInformation="none")
  if(!is.null(prepared$categorical)) result$contrasts <- prepared$categorical$contrasts
  if(!is.null(curveRandom)){
    result$curveRandomMode <- stats::setNames(
      result$randomMode[curveRandom$start+seq_len(curveRandom$count)],curveRandom$columns)
    result$additiveRandomMode <- result$randomMode[seq_len(curveRandom$start)]
  }
  class(result) <- c("mmes.nonlinear","mmes")
  .nonlinear_fit_diagnostics(result,optimized,control,profileStart,preparationSeconds)
}

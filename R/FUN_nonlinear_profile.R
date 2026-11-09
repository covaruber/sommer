.nonlinear_fit_diagnostics <- function(object, optimized, control,
    started=NULL, preparation=NULL){
  diagnostics <- object$nonlinearDiagnostics
  diagnostics$optimizer <- list(method="nlminb",gradient="analytic",
    iterations=object$optimizer$iterations,evaluations=object$optimizer$evaluations,
    message=object$optimizer$message)
  diagnostics$starts <- list(curve=if(is.null(control$task$nonlinear$start)) "auto" else "supplied",
    initialCoordinates=diagnostics$initialParameters,
    conditional=if(!isTRUE(object$nonlinearFamily) && is.null(object$nonlinearRandom))
      "exact Gaussian solve" else "zero and +/-0.5 deterministic latent coordinates")
  diagnostics$dimensions <- c(observations=object$nobs,predictors=length(object$input),
    fixed=ncol(object$fixedDesign),curve=length(object$nonlinearCoordinates),
    randomFull=length(object$randomMode),randomLatent=length(object$latentRandomMode),
    freeParameters=length(object$optimizer$par),
    fixedNonzeros=Matrix::nnzero(object$fixedDesign),
    randomNonzeros=Matrix::nnzero(object$randomDesign))
  diagnostics$memory <- list(budgetMB=control$execution$denseMemoryMB,
    sparseDesignBytes=as.numeric(utils::object.size(object$fixedDesign))+
      as.numeric(utils::object.size(object$randomDesign)),
    admission="Sparse-support component-wise dense factor-fill bound; checked before Hessian assembly",
    scope="Conservative operation workspace, not a process-wide peak RSS limit")
  if(!is.null(started)){
    total <- proc.time()[["elapsed"]]-started
    diagnostics$profile <- list(seconds=c(preparation=preparation,
      outerOptimization=unname(optimized$profile["outerTotal"]),
      resultConstruction=max(0,total-preparation-optimized$profile["outerTotal"])),
      marginalEvaluationSeconds=unname(optimized$profile["marginalEvaluations"]),
      totalSeconds=total,
      note="Marginal evaluation time is included in outerOptimization, not additive.")
  }
  object$nonlinearDiagnostics <- diagnostics
  object
}

profile_nonlinear <- function(object, replications=3L, control=NULL,
    inference=FALSE){
  if(!inherits(object,"mmes.nonlinear"))
    stop("profile_nonlinear() requires a nonlinear mmes fit.",call.=FALSE)
  replications <- .mme_positive(replications,"replications",integer=TRUE)
  .mme_flag(inference,"inference")
  settings <- .nonlinear_inference_settings(object,control)
  object$control$execution$denseMemoryMB <- settings$memory
  scope <- .mmes_blas_scope(object$control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  creation <- system.time(context <- .nonlinear_postfit_context(object,settings$memory))[["elapsed"]]
  parameters <- object$optimizer$par
  value <- context$evaluate(parameters,gradient=FALSE)
  if(abs(value$logLik-object$llik)>1e-7*(1+abs(object$llik)))
    stop("Profile likelihood differs from the saved fit.",call.=FALSE)
  curve <- object$nonlinear
  operations <- list(
    curveMean=function() sommer_nonlinear_curve_eval(curve$curve,object$input,
      object$nonlinearCoordinates,gradient=FALSE),
    curveJacobian=function() sommer_nonlinear_curve_eval(curve$curve,object$input,
      object$nonlinearCoordinates,gradient=TRUE),
    conditionalCurvature=function() context$conditional(parameters,as.numeric(value$latentMode)),
    marginalObjective=function() context$evaluate(parameters,gradient=FALSE),
    marginalScore=function() context$evaluate(parameters,gradient=TRUE),
    conditionalPrediction=function() stats::predict(object))
  if(inference){
    local <- object
    local$nonlinearInference <- NULL
    operations$postFitInference <- function() post_nonlinear(local,control=control)
  }
  table <- do.call(rbind,lapply(names(operations),function(phase){
    operation <- operations[[phase]]
    operation()
    elapsed <- numeric(replications)
    for(index in seq_len(replications)){
      start <- proc.time()[["elapsed"]]
      operation()
      elapsed[index] <- proc.time()[["elapsed"]]-start
    }
    data.frame(phase=phase,seconds=stats::median(elapsed),
      minimumSeconds=min(elapsed),replications=replications)
  }))
  diagnostics <- object$nonlinearDiagnostics
  list(table=table,contextSeconds=unname(creation),
    fitProfile=diagnostics$profile,dimensions=diagnostics$dimensions,
    work=context$evaluate(parameters)$work,
    optimizer=diagnostics$optimizer,starts=diagnostics$starts,
    memory=c(diagnostics$memory,list(fitBytes=as.numeric(utils::object.size(object)),
      profileBudgetMB=settings$memory)),
    likelihood=object$likelihoodMethod,converged=isTRUE(object$convergence),
    boundaryParameters=diagnostics$boundaryParameters,
    environment=list(R=R.version.string,platform=R.version$platform,
      BLAS=unname(extSoftVersion()["BLAS"]),
      threads=object$control$execution$blasThreads),
    note=paste("Warmup excluded; operation times are inclusive, overlapping and not additive.",
      "conditionalCurvature fixes the latent mode; marginalObjective includes mode search;",
      "marginalScore additionally includes analytic derivative contractions.",
      "Memory reports retained R object sizes and admission policy, not native peak RSS."))
}

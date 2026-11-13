.spec_profile_fit <- function(object){
  if(!inherits(object,"mmes.specification") || !is.null(object$baseFit))
    stop("This method requires a reference mmes.specification fit.",call.=FALSE)
  if(!isTRUE(object$convergence))
    stop("Specification profiling requires a converged fit.",call.=FALSE)
  if(!object$likelihoodMethod %in% c("exact ML","Laplace ML"))
    stop("Specification profiling requires an exact ML or Laplace ML fit.",call.=FALSE)
  invisible(object)
}

.spec_profile <- function(object,replications=3L,control=NULL,inference=FALSE){
  if(!inherits(object,"mmes.specification") ||
     is.null(object$likelihoodContext) || is.null(object$optimizer$par))
    stop("profile_nonlinear() requires a reference mmes.specification fit.",
      call.=FALSE)
  replications <- .mme_positive(replications,"replications",integer=TRUE)
  .mme_flag(inference,"inference")
  settings <- .nonlinear_inference_settings(object,control)
  scope <- .mmes_blas_scope(object$control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  context <- object$likelihoodContext
  parameters <- object$optimizer$par
  if(context$p)
    .nonlinear_workspace(as.double(context$n)*context$p,settings$memory,
      "Specification profile Jacobian")
  if(context$q)
    .nonlinear_workspace(as.double(context$q)^2,settings$memory,
      "Specification profile conditional curvature")
  creation <- 0
  value <- .spec_evaluate(parameters,context,details=TRUE)
  if(is.null(value) || !is.finite(value$value) ||
     abs(-value$value-object$llik)>1e-7*(1+abs(object$llik)))
    stop("Profile objective differs from the saved specification fit.",call.=FALSE)
  mode <- unlist(value$mode,use.names=FALSE)
  splitMode <- .spec_split_mode(mode,context$widths)
  predictor <- function() .spec_predictor(context$spec,
    parameters[seq_len(context$p)],context$data,context$design,
    context$randomDesign,splitMode)
  operations <- list(predictor=predictor)
  if(context$p){
    predictorFromCoefficients <- function(coefficients)
      .spec_predictor(context$spec,coefficients,context$data,context$design,
        context$randomDesign,splitMode)
    operations$jacobian <- function()
      .nonlinear_jacobian(predictorFromCoefficients,
        parameters[seq_len(context$p)],settings$numerical,
        "Specification predictor Jacobian")
  }
  if(context$q){
    conditional <- .spec_conditional(parameters,context,value$covariances)
    operations$conditionalCurvature <- function()
      .nonlinear_jacobian(conditional$gradient,mode,settings$numerical,
        "Specification conditional Hessian")
  }
  objective <- function(at) .spec_evaluate(at,context)
  operations$objective <- function() objective(parameters)
  operations$score <- function() .spec_gradient(parameters,objective)
  operations$prediction <- function() stats::predict(object)
  if(inference){
    local <- object
    local$nonlinearInference <- NULL
    operations$postFitInference <- function()
      .spec_post_nonlinear(local,mode=1L,control=control)
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
  optimizer <- object$optimizer
  memory <- list(budgetMB=settings$memory,
    retainedDesignBytes=as.numeric(utils::object.size(context$design))+
      as.numeric(utils::object.size(context$randomDesign)),
    fitBytes=as.numeric(utils::object.size(object)),
    profileBudgetMB=settings$memory,
    admission="Reference dense workspaces are admitted during model preparation.")
  list(table=table,contextSeconds=creation,
    fitProfile=object$specificationDiagnostics$profile,
    dimensions=c(observations=context$n,meanCoefficients=context$p,
      randomCoordinates=context$q,randomRecords=length(context$records),
      freeParameters=length(parameters)),
    work=list(outerScore=object$specificationDiagnostics$outerScore,
      modeScore=object$modeScore,competingModes=object$competingModes,
      failedModeStarts=object$failedModeStarts),
    optimizer=list(method="nlminb",
      gradient="numerical finite differences",
      iterations=optimizer$iterations,evaluations=optimizer$evaluations,
      message=optimizer$message),
    starts=list(method="deterministic reference initialization; no refit"),
    memory=memory,likelihood=object$likelihoodMethod,
    converged=isTRUE(object$convergence),
    boundaryParameters=object$specificationDiagnostics$boundaryParameters,
    environment=list(R=R.version.string,platform=R.version$platform,
      BLAS=unname(extSoftVersion()["BLAS"]),
      threads=object$control$execution$blasThreads),
    note=paste("Warmup excluded; operation times are inclusive, overlapping and not additive.",
      "The saved likelihood context and fitted coordinates are reused; the outer model is not refit.",
      "Objective and score operations may optimize conditional modes as required by the likelihood.",
      "The reference evaluator reports no native cumulative-work counters."))
}

.spec_comparison_contract <- function(object){
  context <- object$likelihoodContext
  parameters <- lapply(context$spec$parameters,function(parameter){
    parameter$positions <- NULL
    parameter
  })
  list(data=context$data,response=context$y,likelihoodMethod=object$likelihoodMethod,
    familySpecifications=context$familySpecifications,
    familyIndex=context$familyIndex,trials=context$trials,
    parameters=parameters,
    expressions=context$spec$expressions,inputs=context$spec$inputs,
    parameterDesign=context$design$parameters,records=context$records,
    randomDesign=context$randomDesign,exact=context$exact)
}

.spec_nested_fixed_terms <- function(smaller,larger){
  smallTerms <- stats::terms(smaller$fixedFormula)
  largeTerms <- stats::terms(larger$fixedFormula)
  smallLabels <- attr(smallTerms,"term.labels")
  largeLabels <- attr(largeTerms,"term.labels")
  if(attr(smallTerms,"intercept")!=attr(largeTerms,"intercept") ||
     !all(smallLabels %in% largeLabels) ||
     identical(smallLabels,largeLabels))
    return(FALSE)
  smallDesign <- smaller$likelihoodContext$design$fixed
  largeDesign <- larger$likelihoodContext$design$fixed
  common <- colnames(smallDesign$x)
  if(!all(common %in% colnames(largeDesign$x)) ||
     !identical(smallDesign$offset,largeDesign$offset) ||
     !isTRUE(all.equal(smallDesign$x,
       largeDesign$x[,match(common,colnames(largeDesign$x)),drop=FALSE],
       check.attributes=FALSE)))
    return(FALSE)
  added <- setdiff(colnames(largeDesign$x),common)
  length(larger$optimizer$par)-length(smaller$optimizer$par)==length(added)
}

.spec_wald <- function(object,parm,null){
  .spec_profile_fit(object)
  if(missing(parm) || missing(null) || !is.character(parm) ||
     !length(parm) || anyNA(parm) || any(!nzchar(parm)) ||
     anyDuplicated(parm) || any(!parm %in% names(object$b)) ||
     !is.numeric(null) || is.complex(null) || length(null)!=length(parm) ||
     any(!is.finite(null)) || is.null(names(null)) || anyNA(names(null)) ||
     anyDuplicated(names(null)) || !setequal(names(null),parm))
    stop("A Wald test requires explicit mean-coordinate parm names and one finite, exactly named null per coordinate.",
      call.=FALSE)
  null <- null[parm]
  object <- .spec_post_nonlinear(object,0L)
  covariance <- stats::vcov(object,type="mean")[parm,parm,drop=FALSE]
  difference <- object$b[parm]-null
  statistic <- as.numeric(crossprod(difference,solve(covariance,difference)))
  result <- data.frame(Df=length(parm),Chisq=statistic,
    `Pr(>Chisq)`=stats::pchisq(statistic,length(parm),lower.tail=FALSE),
    check.names=FALSE)
  class(result) <- c("anova","data.frame")
  attr(result,"heading") <- paste("Wald test for explicit finite null values on regular mean coordinates:",
    paste(parm,collapse=", "),".")
  result
}

.spec_likelihood_compare <- function(smaller,larger){
  .spec_profile_fit(smaller)
  .spec_profile_fit(larger)
  if(!identical(.spec_comparison_contract(smaller),
      .spec_comparison_contract(larger)))
    stop("Specification likelihood comparisons require identical data, nonlinear expressions, parameter designs, random covariance structures, families, and likelihood methods.",
      call.=FALSE)
  if(!.spec_nested_fixed_terms(smaller,larger))
    stop("Specification comparisons require strictly nested ordinary fixed terms with unchanged nonlinear and covariance structure; supply the smaller model first.",
      call.=FALSE)
  .spec_post_nonlinear(smaller,0L)
  .spec_post_nonlinear(larger,0L)
  degrees <- length(larger$optimizer$par)-length(smaller$optimizer$par)
  statistic <- 2*(larger$llik-smaller$llik)
  tolerance <- 1e-7*(1+abs(smaller$llik)+abs(larger$llik))
  if(degrees<=0L || statistic < -tolerance)
    stop("The larger fixed-effects model must have a non-decreasing fitted likelihood and positive added degrees of freedom.",
      call.=FALSE)
  statistic <- max(0,statistic)
  result <- data.frame(npar=c(length(smaller$optimizer$par),
      length(larger$optimizer$par)),
    logLik=c(smaller$llik,larger$llik),
    Chisq=c(NA_real_,statistic),Df=c(NA_real_,degrees),
    `Pr(>Chisq)`=c(NA_real_,stats::pchisq(statistic,degrees,lower.tail=FALSE)),
    check.names=FALSE)
  rownames(result) <- c("smaller","larger")
  class(result) <- c("anova","data.frame")
  attr(result,"heading") <- paste("Likelihood-ratio test for nested ordinary fixed terms;",
    "all nonlinear, parameter, data, covariance and family semantics are identical.")
  result
}

anova.mmes.specification <- function(object,object2=NULL,...,parm=NULL,null=NULL){
  .spec_dots(...)
  if(is.null(object2)){
    if(is.null(parm) || is.null(null))
      stop("Supply object2 for a nested fixed-term comparison, or explicit parm and null for a mean-coordinate Wald test.",
        call.=FALSE)
    return(.spec_wald(object,parm,null))
  }
  if(!is.null(parm) || !is.null(null))
    stop("Specify either object2 or a Wald-test parm/null, not both.",call.=FALSE)
  .spec_likelihood_compare(object,object2)
}

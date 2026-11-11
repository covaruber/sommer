#' Model fitting method
#'
#' Declare the statistical objective independently of numerical controls.
#' Defaults preserve Gaussian REML and non-Gaussian mixed-model PQL with
#' REML Gaussian working fits. This is not an exact GLMM REML likelihood.
#' Fixed-only generalized and nonlinear likelihood models require an explicit
#' \code{objective="ml"}.
#' @param objective \code{"reml"} or \code{"ml"}.
#' @param integration \code{"auto"}, \code{"exact"}, or \code{"laplace"}.
#' @param approximation \code{"auto"}, \code{"pql"}, or \code{"likelihood"}.
#' @param task \code{"fit"} or the existing known-covariance \code{"solve"} task.
#' @export
mmes_method <- function(objective="reml", integration="auto",
    approximation="auto", task="fit"){
  structure(list(
    objective=.mme_choice(objective,c("reml","ml"),"method$objective"),
    integration=.mme_choice(integration,c("auto","exact","laplace"),
      "method$integration"),
    approximation=.mme_choice(approximation,c("auto","pql","likelihood"),
      "method$approximation"),
    task=.mme_choice(task,c("fit","solve"),"method$task")),
    class="sommer_mmes_method")
}

#' Modular numerical controls
#'
#' Statistical objectives belong in \code{mmes_method()}, parameter starts in
#' \code{mmes(start=)}, and family parameter choices in \code{mmes(familyParameters=)}.
#' Existing equation, solver, covariance-update, execution and output modules
#' retain their mathematical responsibilities.
#' @param optimizer Outer optimizer settings: \code{maxit}, \code{gradient},
#'   and \code{algorithm}. Only \code{"nlminb"} is currently supported for
#'   nonlinear/Laplace optimization; Gaussian covariance engines retain their
#'   existing AI/update modules.
#' @param mode Conditional mode settings: \code{maxit}, \code{tol},
#'   \code{stepTol}, \code{warmStart}, and \code{starts}.
#' @param working PQL working-fit iteration settings.
#' @param likelihood Native likelihood algebra settings: \code{backend},
#'   \code{trace}, and \code{blocks}.
#' @param population Population-prediction quadrature settings.
#' @param inference Post-fit nonlinear inference settings.
#' @param formulation,representation,solve,logdet,score,curvature,update,proposal,convergence,output,execution,diagnostics
#'   Existing modular numerical settings; see the numerical controls vignette.
#' @export
mmes_control <- function(formulation="henderson", optimizer=list(), mode=list(),
    working=list(), likelihood=list(), population=list(), inference=list(),
    representation=list(), solve=list(), logdet=list(), score=list(),
    curvature=list(), update=list(), proposal=list(), convergence=list(),
    output=list(), execution=list(), diagnostics=list()){
  optimizer <- .mme_module(optimizer,
    list(maxit=NULL,gradient="auto",algorithm="nlminb"),"optimizer")
  if(!is.null(optimizer$maxit))
    optimizer$maxit <- .mme_positive(optimizer$maxit,"optimizer$maxit",TRUE)
  optimizer$gradient <- .mme_choice(optimizer$gradient,
    c("auto","analytic","central"),"optimizer$gradient")
  optimizer$algorithm <- .mme_choice(optimizer$algorithm,"nlminb",
    "optimizer$algorithm")
  mode <- .mme_module(mode,list(maxit=100L,tol=1e-8,stepTol=1e-8,
    warmStart="auto",starts="deterministic"),"mode")
  mode$maxit <- .mme_positive(mode$maxit,"mode$maxit",TRUE)
  for(name in c("tol","stepTol"))
    mode[[name]] <- .mme_positive(mode[[name]],paste0("mode$",name))
  mode$warmStart <- .mme_choice(mode$warmStart,c("off","auto","on"),"mode$warmStart")
  mode$starts <- .mme_choice(mode$starts,"deterministic","mode$starts")
  working <- .mme_module(working,.mme_pql_defaults,"working")
  likelihood <- .mme_module(likelihood,list(backend="auto",
    trace="panels",blocks="off"),"likelihood")
  population <- .mme_module(population,.mme_nonlinear_prediction_defaults,
    "population")
  inference <- .mme_module(inference,.mme_nonlinear_inference_defaults,"inference")
  checked <- .mmes_control_build(task=list(pql=working,glmm=likelihood,
    nonlinear=list(population=population,inference=inference)))
  working <- checked$task$pql
  population <- checked$task$nonlinear$population
  inference <- checked$task$nonlinear$inference
  args <- list(formulation=formulation,representation=representation,solve=solve,
    logdet=logdet,score=score,curvature=curvature,update=update,proposal=proposal,
    convergence=convergence,output=output,execution=execution,
    diagnostics=diagnostics)
  internal <- do.call(.mmes_control_internal,args)
  settings <- list(optimizer=optimizer,mode=mode,
    working=working,likelihood=likelihood,population=population,inference=inference)
  result <- c(unclass(internal)[names(args)],settings)
  class(result) <- "sommer_mmes_control"
  attr(result,"modelControls") <- settings
  supplied <- as.list(match.call())[-1L]
  attr(result,"publicSupplied") <- names(supplied)
  attr(result,"publicArguments") <- unclass(result)
  attr(result,"explicitMaxIters") <- "maxIters" %in% names(convergence)
  result
}

.mmes_engine_arguments <- function(control){
  value <- unclass(control)
  value[intersect(names(value),names(formals(.mmes_control_build)))]
}

.mmes_resolve_controls <- function(control,method,nonlinear,start,familyParameters){
  if(!inherits(control,"sommer_mmes_control") ||
     is.null(attr(control,"modelControls")))
    stop("control must be created by mmes_control().",call.=FALSE)
  if(!inherits(method,"sommer_mmes_method") ||
     !identical(unclass(method),unclass(do.call(mmes_method,unclass(method)))))
    stop("method must be created by mmes_method().",call.=FALSE)
  supplied <- attr(control,"publicSupplied")
  explicitMaxIters <- isTRUE(attr(control,"explicitMaxIters")) ||
    control$convergence$maxIters!=.mmes_control_build()$convergence$maxIters
  checked <- do.call(mmes_control,unclass(control))
  attr(checked,"publicSupplied") <- supplied
  control <- checked
  settings <- attr(control,"modelControls")
  args <- .mmes_engine_arguments(control)
  args$task <- list(mode=method$task,objective=method$objective,
    pql=settings$working,
    glmm=c(list(method=if(!nonlinear && method$integration=="laplace") "laplace" else "pql",
      familyParameters=familyParameters),settings$likelihood))
  args$task$glmm$warmStart <- if(nonlinear) "auto" else settings$mode$warmStart
  if(!is.null(settings$optimizer$maxit)){
    if(explicitMaxIters &&
       control$convergence$maxIters != settings$optimizer$maxit)
      stop("Specify the outer iteration limit once: optimizer$maxit or convergence$maxIters.",
        call.=FALSE)
    args$convergence$maxIters <- settings$optimizer$maxit
    if(!nonlinear) args$task$glmm$maxit <- settings$optimizer$maxit
  }
  if(!nonlinear){
    args$task$glmm$modeMaxit <- settings$mode$maxit
    args$task$glmm$modeTol <- settings$mode$tol
    args$task$glmm$gradient <- if(settings$optimizer$gradient=="central")
      "numeric" else "analytic"
  }
  if(nonlinear){
    if(settings$mode$warmStart=="on")
      stop("Nonlinear conditional modes currently require mode$warmStart='auto' or 'off'; auto selects deterministic cold starts.",
        call.=FALSE)
    args$task$nonlinear <- list(
      start=if(is.list(start) && "parameters" %in% names(start))
        start$parameters else start,
      mode=list(maxIters=settings$mode$maxit,scoreTol=settings$mode$tol,
        stepTol=settings$mode$stepTol,warmStart="off",starts=settings$mode$starts),
      population=settings$population,inference=settings$inference)
  }else if(!is.null(start)){
    stop("start is currently supported for declared nonlinear predictors only; use covPar or covariance descriptors for covariance initialization.",
      call.=FALSE)
  }
  if(method$approximation=="pql" && method$integration!="auto")
    stop("PQL is a working-model approximation, not exact or Laplace likelihood integration.",
      call.=FALSE)
  if(!nonlinear && method$approximation=="likelihood" && method$integration=="auto")
    args$task$glmm$method <- "laplace"
  if(method$approximation=="pql" && nonlinear)
    stop("PQL is not supported for nonlinear parameter models.",call.=FALSE)
  result <- do.call(.mmes_control_internal,args)
  attr(result,"publicControl") <- control
  result
}

#' Unified model specification
#'
#' The formula describes a link-scale predictor. Named \code{parameter()}
#' declarations distinguish unknown nonlinear parameters from data columns.
#' Named \code{predictors} can contain native curves or bounded graphs and are
#' referenced with \code{graph_term("name")} in the formula.
#' @param fixed,random,rcov,data,W,weights,naMethodX,naMethodY,naMethodRandom,naMethodR,contrasts,family,vcc,covPar
#'   Statistical model and observation-preparation inputs.
#' @param parameters Named nonlinear parameter declarations.
#' @param predictors Named native curve or graph components.
#' @param start Natural-scale nonlinear parameter starts or a list containing
#'   \code{parameters}; coefficient models use explicitly named coefficient starts.
#' @param method A validated \code{mmes_method()}.
#' @param control A validated \code{mmes_control()}.
#' @param familyParameters Fixed/estimated family nuisance parameter declarations.
#' @export
mmes <- function(fixed,random=NULL,rcov=NULL,data,W,weights=NULL,
    naMethodX="exclude",naMethodY="exclude",naMethodRandom="exclude",
    naMethodR="exclude",contrasts=NULL,family=stats::gaussian(),vcc=NULL,
    covPar=NULL,parameters=list(),predictors=list(),start=NULL,
    method=mmes_method(),control=mmes_control(),familyParameters=list()){
  if(!inherits(fixed,"formula") || length(fixed)!=3L)
    stop("fixed must be a two-sided model formula.",call.=FALSE)
  if(!is.data.frame(data)) stop("data must be a data frame.",call.=FALSE)
  if(!is.list(parameters) || !is.list(predictors))
    stop("parameters and predictors must be named lists.",call.=FALSE)
  for(name in c("parameters","predictors")){
    value <- get(name)
    if(length(value) && (is.null(names(value)) || anyNA(names(value)) ||
        any(!nzchar(names(value))) || anyDuplicated(names(value))))
      stop(name," must have unique non-empty names.",call.=FALSE)
  }
  nonlinear <- length(parameters)>0L || length(predictors)>0L
  if(is.list(start) && (is.null(names(start)) || anyNA(names(start)) ||
      anyDuplicated(names(start)) ||
      any(!names(start) %in% c("parameters","coefficients"))))
    stop("start must contain uniquely named parameters and/or coefficients.",
      call.=FALSE)
  internal <- .mmes_resolve_controls(control,method,nonlinear,start,familyParameters)
  if(nonlinear && (method$objective!="ml" || method$task!="fit"))
    stop("Nonlinear predictors require method=mmes_method(objective='ml', task='fit').",
      call.=FALSE)
  if(!nonlinear && internal$task$glmm$method=="laplace" &&
     method$objective!="ml")
    stop("Laplace likelihood fitting requires method$objective='ml'; REML remains available for Gaussian and PQL working fits.",
      call.=FALSE)
  gaussian <- inherits(family,"family") && identical(family$family,"gaussian") &&
    identical(family$link,"identity")
  if(!nonlinear && is.null(random) && !gaussian && method$objective!="ml")
    stop("Fixed-only generalized likelihood models require method$objective='ml'; REML/PQL defaults apply to mixed working models.",
      call.=FALSE)
  if(!nonlinear && !gaussian && !is.null(random) &&
     method$integration=="exact")
    stop("Exact integration is not available for non-Gaussian mixed models.",
      call.=FALSE)
  if(!nonlinear && !is.null(random) && internal$task$glmm$method=="pql" &&
     !identical(control$mode,
       attr(mmes_control(),"modelControls")$mode))
    stop("mode controls conditional Laplace modes, not Gaussian or PQL working fits.",
      call.=FALSE)
  if(!nonlinear && internal$task$glmm$method=="laplace" &&
     control$mode$stepTol!=1e-8)
    stop("mode$stepTol is currently supported by nonlinear modes only; native GLMM modes use mode$tol.",
      call.=FALSE)
  args <- list(fixed=fixed,data=data,weights=weights,
    naMethodX=naMethodX,naMethodY=naMethodY,naMethodRandom=naMethodRandom,
    naMethodR=naMethodR,contrasts=contrasts,family=family,vcc=vcc,
    covPar=covPar,control=internal)
  if(!is.null(random)) args$random <- random
  if(!missing(W)) args$W <- W
  if(!is.null(rcov)) args$rcov <- rcov
  native <- .mmes_native_predictor(fixed,parameters,predictors)
  if(!is.null(native) && is.null(random) &&
     !length(native$model$random) && is.null(native$model$joint)){
    if(inherits(native$model,"sommer_nonlinear_curve")){
      declarations <- lapply(seq_along(native$model$parameters),function(i)
        parameter(link=if(native$model$positive[i]) "log" else "identity"))
      names(declarations) <- native$model$parameters
      parameters <- declarations
      args$fixed[[3L]] <- call("+",native$fixed[[3L]],
        .mmes_curve_expression(native$model))
      values <- if(is.list(start)) start$parameters else start
      qualified <- paste0(native$name,"::",native$model$parameters)
      if(!is.null(values) && setequal(names(values),qualified)){
        values <- values[qualified]
        names(values) <- native$model$parameters
        start <- list(parameters=values)
      }
      if(is.null(start)){
        frame <- stats::model.frame(native$fixed,data,
          na.action=stats::na.fail)
        response <- stats::model.response(frame)
        offset <- stats::model.offset(frame)
        if(is.null(offset)) offset <- 0
        theta <- .nonlinear_start(native$model,data[[native$model$input]],
          if(gaussian) response-offset else rep(1,nrow(data)),
          internal$task$nonlinear)
        start <- list(parameters=.nonlinear_natural(native$model,theta))
      }
      predictors <- list()
    }else{
      values <- if(is.list(start)) start$parameters else start
      if(!is.null(values) && setequal(names(values),native$model$parameters)){
        names(values) <- paste0(native$name,"::",names(values))
        start <- list(parameters=values)
      }
    }
    native <- NULL
  }
  if(!is.null(native)){
    if(is.list(start) && "coefficients" %in% names(start))
      stop("Native graph/curve components require start$parameters, not coefficient starts.",
        call.=FALSE)
    args$fixed <- native$fixed
    args$nonlinear <- native$model
    requestedGradient <- control$optimizer$gradient
    graphPredictor <- inherits(native$model,"sommer_nonlinear_graph")
    numericalGraph <- graphPredictor &&
      (!.nonlinear_exact_gaussian(family,internal) ||
       length(native$model$random)>0L || !is.null(native$model$joint))
    if((requestedGradient=="analytic" && numericalGraph) ||
       (requestedGradient=="central" && !numericalGraph))
      stop("Requested optimizer gradient is unavailable for this native predictor; use gradient='auto' or the supported central/analytic method.",
        call.=FALSE)
    if(!is.null(internal$task$nonlinear$start)){
      coordinates <- internal$task$nonlinear$start
      qualified <- paste0(native$name,"::",native$model$parameters)
      if(setequal(names(coordinates),qualified)){
        coordinates <- coordinates[qualified]
        names(coordinates) <- native$model$parameters
        internal$task$nonlinear$start <- coordinates
        args$control <- internal
      }
    }
    if(method$integration=="exact" &&
       (!gaussian || length(native$model$random) || !is.null(native$model$joint)))
      stop("Exact integration is unavailable for this nonlinear random-effect model.",
        call.=FALSE)
    if(inherits(native$model,"sommer_nonlinear_graph") &&
       !is.null(internal$task$nonlinear$start)){
      coordinates <- internal$task$nonlinear$start
      positive <- vapply(native$model$parameterMap,function(x)
        identical(x$transform,"log"),logical(1))
      if(any(!native$model$parameters %in% names(coordinates)))
        stop("start must name every graph coordinate.",call.=FALSE)
      if(any(coordinates[native$model$parameters[positive]]<=0))
        stop("Log-linked graph coordinates require positive natural-scale starts.",
          call.=FALSE)
      coordinates[native$model$parameters[positive]] <-
        log(coordinates[native$model$parameters[positive]])
      internal$task$nonlinear$start <- coordinates
      args$control <- internal
    }
    result <- do.call(.mmes_dispatch,args)
  }else if(nonlinear || (is.null(random) && !gaussian)){
    args$parameters <- parameters
    args$predictors <- predictors
    args$start <- start
    args$specificationMethod <- unclass(method)
    args$familyParameters <- familyParameters
    result <- do.call(.mmes_specification_fit,args)
  }else{
    result <- do.call(.mmes_dispatch,args)
  }
  result$method <- method
  result$.numericalControl <- result$control
  result$control <- control
  if(is.null(result$executionPlan))
    result$executionPlan <- list(version=2L,task=method$task,
      objective=result$likelihoodMethod,derivatives=result$derivativeMethod,
      integration=if(grepl("Laplace",result$likelihoodMethod,fixed=TRUE))
        "laplace" else "exact",
      status=if(isTRUE(result$convergence)) "converged" else "not converged")
  result$executionPlan$requestedMethod <- unclass(method)
  result$call <- match.call()
  result
}

.mmes_native_predictor <- function(fixed,parameters,predictors){
  if(length(parameters) || !length(predictors)) return(NULL)
  if(length(predictors)!=1L) return(NULL)
  model <- predictors[[1L]]
  if(!inherits(model,"sommer_nonlinear_graph") &&
     !inherits(model,"sommer_nonlinear_curve"))
    stop("predictors must contain nonlinear_graph() or nonlinear_curve() specifications.",
      call.=FALSE)
  terms <- .mmes_split_plus(fixed[[3L]])
  marker <- vapply(terms,function(term) is.call(term) &&
    identical(term[[1L]],as.name("graph_term")),logical(1))
  if(sum(marker)!=1L || length(terms[[which(marker)]])!=2L ||
     !identical(terms[[which(marker)]][[2L]],names(predictors)[1L]))
    return(NULL)
  ordinary <- terms[!marker]
  fixed[[3L]] <- if(!length(ordinary)) 1 else
    Reduce(function(a,b) call("+",a,b),ordinary)
  list(fixed=fixed,model=model,name=names(predictors)[1L])
}

.mmes_prediction_arguments <- function(arguments,target){
  if(is.null(target)) return(arguments)
  target <- .mme_choice(target,c("conditional","zero","population"),"target")
  if("re.form" %in% names(arguments))
    stop("Specify target or re.form, not both.",call.=FALSE)
  arguments["re.form"] <- list(switch(target,conditional=NULL,zero=NA,
    population="population"))
  arguments
}

.mmes_curve_expression <- function(model){
  input <- as.name(model$input)
  switch(model$curve,
    exponential=substitute(asymptote*exp(-rate*x),list(x=input)),
    logistic=substitute(asymptote/(1+exp(-rate*(x-midpoint))),list(x=input)),
    saturation=substitute(asymptote*x/(half+x),list(x=input)))
}

.mmes_fitting_control <- function(object){
  control <- object$.numericalControl
  if(is.null(control)) return(object$control)
  control$execution <- object$control$execution
  control$output <- object$control$output
  control$task$nonlinear$population <- object$control$population
  control$task$nonlinear$inference <- object$control$inference
  control
}

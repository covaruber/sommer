.nonlinear_check_dots <- function(...){
  if(length(list(...))){
    stop("Unused arguments are not supported for nonlinear models.", call.=FALSE)
  }
}

coef.mmes.nonlinear <- function(object, ..., type=c("mean","working","all")){
  .nonlinear_check_dots(...)
  type <- match.arg(type)
  if(type=="working") return(object$optimizer$par)
  if(type=="all") return(.nonlinear_parameter_values(object,object$optimizer$par))
  c(if(length(object$b)) stats::setNames(object$b, paste0("fixed:", names(object$b))),
    stats::setNames(object$nonlinearParameters,
      paste0("curve:", names(object$nonlinearParameters))))
}

fitted.mmes.nonlinear <- function(object, type=c("response", "link"), ...){
  .nonlinear_check_dots(...)
  type <- match.arg(type)
  if(isTRUE(object$nonlinearFamily))
    return(.nonlinear_family_output(object,object$linear.predictors,type))
  object$fitted.values
}

residuals.mmes.nonlinear <- function(object, type="response", ...){
  .nonlinear_check_dots(...)
  if(isTRUE(object$nonlinearFamily)){
    proxy <- object
    class(proxy) <- c("mmes.glmm.laplace","mmes.glmm","mmes")
    value <- residuals(proxy,type=type)
    return(if(object$singleCategorical) value[[1L]] else value)
  }
  if(!identical(type, "response")){
    stop("Nonlinear models currently support only response residuals.", call.=FALSE)
  }
  object$y-object$fitted.values
}

logLik.mmes.nonlinear <- function(object, ...){
  .nonlinear_check_dots(...)
  structure(object$llik, class="logLik", nobs=object$nobs, df=object$npar)
}

nobs.mmes.nonlinear <- function(object, ...){
  .nonlinear_check_dots(...)
  object$nobs
}

summary.mmes.nonlinear <- function(object, ..., inference=FALSE){
  .nonlinear_check_dots(...)
  .mme_flag(inference,"inference")
  uncertainty <- if(inference) .nonlinear_uncertainty(object) else object$nonlinearInference
  structure(list(curve=object$nonlinear$curve, coefficients=stats::coef(object),
    parameterTable=if(!is.null(uncertainty)) data.frame(
      Estimate=uncertainty$parameters,
      StdError=sqrt(pmax(0,diag(uncertainty$parameterCovariance)))) else NULL,
    uncertaintyMethod=uncertainty$method,
    groupParameterSE=uncertainty$groupParameterSE,
    covariance=object$covPar,
    residualVariance=if(!isTRUE(object$nonlinearFamily)) unname(object$dispersion) else NULL,
    families=if(isTRUE(object$nonlinearFamily)) vapply(object$familySpecifications,
      function(specification) paste0(specification$family,"(",specification$link,")"),
      character(1)) else NULL,
    familyParameters=if(isTRUE(object$nonlinearFamily))
      lapply(object$familySpecifications,function(specification){
        if(specification$family=="ordinal") return(specification$thresholds)
        if(!is.null(specification$zeroDesign)) return(specification$zeroCoefficients)
        if(!is.null(specification$parameter))
          stats::setNames(.laplace_family_parameter_value(specification),specification$parameter)
      }) else NULL,
    logLik=object$llik, convergence=object$convergence,
    likelihoodMethod=if(is.null(object$likelihoodMethod)) "exact Gaussian ML" else object$likelihoodMethod,
    groupParameters=object$groupParameters,
    diagnostics=object$nonlinearDiagnostics), class="summary.mmes.nonlinear")
}

print.summary.mmes.nonlinear <- function(x, ...){
  cat("Nonlinear mixed model:", x$likelihoodMethod, "\n")
  cat("Curve:", x$curve, "| Converged:", x$convergence, "\n")
  cat("logLik:", format(x$logLik), "\n")
  if(is.null(x$parameterTable)) print(x$coefficients) else print(x$parameterTable)
  if(!is.null(x$residualVariance)) cat("Residual variance:", format(x$residualVariance), "\n")
  if(length(x$families)) cat("Families:", paste(x$families,collapse=", "), "\n")
  if(any(lengths(x$familyParameters))) print(Filter(length,x$familyParameters))
  if(is.null(x$parameterTable)) cat("SE and PEV not computed; request post_nonlinear().\n") else
    cat(x$uncertaintyMethod,"\n")
  invisible(x)
}

print.mmes.nonlinear <- function(x, ...){
  print(summary(x), ...)
  invisible(x)
}

.nonlinear_population_prediction <- function(object, input, design){
  layout <- object$nonlinearRandom
  term <- layout$term
  variance <- .laplace_population_variance(design, object$covStruct[term],
    object$randomRelationships[term], layout$count,
    memoryMB=object$control$execution$denseMemoryMB,
    allowSingular=.nonlinear_has_reduced_rank(object$covStruct[term]))
  settings <- object$control$task$nonlinear$population
  evaluate <- function(order){
    rule <- .laplace_normal_quadrature(order)
    mean <- numeric(length(input))
    for(index in seq_along(rule$nodes)){
      mean <- mean+rule$weights[index]*as.numeric(sommer_nonlinear_curve_random_eval(
        object$nonlinear$curve, input, object$nonlinearCoordinates,
        layout$coordinate, sqrt(variance)*rule$nodes[index]))
    }
    mean
  }
  order <- settings$order
  previous <- evaluate(order)
  repeat {
    order <- min(2L*order, settings$maxOrder)
    current <- evaluate(order)
    if(all(abs(current-previous) <= settings$tol*(1+abs(current)))) return(current)
    if(order >= settings$maxOrder){
      stop("Nonlinear population quadrature did not meet task$nonlinear$population$tol at maxOrder.",
        call.=FALSE)
    }
    previous <- current
  }
}

.nonlinear_random_prediction <- function(object, newdata, design=FALSE){
  variables <- names(object$predictionData)
  if(!all(variables %in% names(newdata))){
    stop("Conditional prediction requires random-effect columns: ",
      paste(setdiff(variables, names(newdata)), collapse=", "), ".", call.=FALSE)
  }
  data <- newdata[variables]
  for(variable in names(object$randomFactorLevels)){
    values <- as.character(data[[variable]])
    levels <- object$randomFactorLevels[[variable]]
    if(anyNA(values) || any(!values %in% levels)){
      stop("Conditional prediction has a missing or unseen random-effect level in '",
        variable, "'.", call.=FALSE)
    }
    data[[variable]] <- factor(values, levels=levels)
  }
  n <- nrow(newdata)
  # Retain the fitted covariance basis even when newdata contains only one group.
  data <- rbind(data, object$predictionData)
  response <- ".sommer_nonlinear_prediction_response"
  while(response %in% names(data)) response <- paste0(response, "_")
  data[[response]] <- 0
  formula <- stats::as.formula(call("~", as.name(response), 0),
    env=environment(object$fixedFormula))
  settings <- object$predictionSettings
  prepared <- .nonlinear_prepare(formula, object$randomFormula, data, object$control,
    "fail", "fail", "fail", "fail", settings$contrasts)
  if(length(prepared$Z) != length(object$randomBlockColumns)){
    stop("Prediction random-effect blocks differ from the fitted model.", call.=FALSE)
  }
  blocks <- lapply(seq_along(prepared$Z), function(index){
    block <- prepared$Z[[index]][seq_len(n), , drop=FALSE]
    fitted <- object$randomBlockColumns[[index]]
    if(is.null(fitted) || is.null(colnames(block)) ||
       anyDuplicated(fitted) || anyDuplicated(colnames(block)) ||
       !setequal(fitted, colnames(block))){
      stop("Prediction random-effect covariance basis differs from the fitted model.",
        call.=FALSE)
    }
    block[, match(fitted, colnames(block)), drop=FALSE]
  })
  combined <- do.call(cbind, blocks)
  if(design) combined else as.numeric(combined %*% object$randomMode)
}

predict.mmes.nonlinear <- function(object, newdata=NULL, re.form=NULL,
    type=c("response", "link"), se.fit=FALSE, cov.fit=FALSE,
    interval=c("none","confidence","prediction"), level=0.95, trials=NULL, ...){
  .nonlinear_check_dots(...)
  type <- match.arg(type)
  interval <- match.arg(interval)
  if(!is.null(trials) && interval!="prediction")
    stop("trials applies only to future-observation prediction intervals.",call.=FALSE)
  if(!identical(se.fit, FALSE) || !identical(cov.fit, FALSE) || interval!="none")
    return(.nonlinear_predict_uncertainty(object,newdata,re.form,type,
      se.fit,cov.fit,interval,level,trials))
  if(isTRUE(object$nonlinearFamily))
    return(.nonlinear_family_predict(object,newdata,re.form,type))
  if(identical(object$nonlinearRandom$version,2L))
    return(.nonlinear_family_predict(.nonlinear_inference_proxy(object),
      newdata,re.form,type))
  includeRandom <- is.null(re.form)
  if(!includeRandom && !identical(re.form, NA) &&
     !identical(re.form, NA_integer_) && !identical(re.form, "population")){
    stop("re.form must be NULL, NA or 'population' for nonlinear models.", call.=FALSE)
  }
  if(is.null(newdata)){
    if(!is.null(object$nonlinearRandom) && identical(re.form, "population")){
      return(object$fixedMean-as.numeric(sommer_nonlinear_curve_eval(
        object$nonlinear$curve, object$input, object$nonlinearCoordinates,
        gradient=FALSE)$mean)+.nonlinear_population_prediction(object, object$input,
          object$nonlinearRandom$design))
    }
    return(if(includeRandom) object$fitted.values else object$fixedMean)
  }
  input <- .nonlinear_input(object$nonlinear, newdata)
  if(!all(object$fixedVariables %in% names(newdata))){
    stop("Prediction requires fixed-effect columns: ",
      paste(setdiff(object$fixedVariables, names(newdata)), collapse=", "), ".",
      call.=FALSE)
  }
  frame <- stats::model.frame(object$fixedTerms, data=newdata,
    xlev=object$fixedXlevels, na.action=stats::na.pass)
  fixed <- stats::model.matrix(object$fixedTerms, frame, contrasts.arg=object$contrasts)
  colnames(fixed)[colnames(fixed) == "(Intercept)"] <- "Intercept"
  if(nrow(fixed) != nrow(newdata) ||
     !identical(colnames(fixed), object$predictionColumns) ||
     any(!is.finite(fixed))){
    stop("Prediction fixed design has missing values or differs from the fitted model.",
      call.=FALSE)
  }
  offset <- stats::model.offset(frame)
  if(is.null(offset)) offset <- numeric(nrow(newdata))
  if(any(!is.finite(offset))){
    stop("Prediction offsets must be finite.", call.=FALSE)
  }
  if(8*as.double(nrow(newdata))*(20+ncol(fixed)) >
     object$control$execution$denseMemoryMB*1048576){
    stop("Nonlinear prediction workspace exceeds execution$denseMemoryMB.", call.=FALSE)
  }
  curveRandom <- object$nonlinearRandom
  if(!is.null(curveRandom) && (includeRandom || identical(re.form, "population"))){
    design <- .nonlinear_random_prediction(object, newdata, design=TRUE)
    curveDesign <- design[, curveRandom$start+seq_len(curveRandom$count), drop=FALSE]
    if(includeRandom){
      deviations <- as.numeric(curveDesign %*% object$curveRandomMode)
      mean <- sommer_nonlinear_curve_random_eval(object$nonlinear$curve, input,
        object$nonlinearCoordinates, curveRandom$coordinate, deviations)
    }else mean <- .nonlinear_population_prediction(object, input, curveDesign)
  }else{
    mean <- sommer_nonlinear_curve_eval(object$nonlinear$curve, input,
      object$nonlinearCoordinates, gradient=FALSE)$mean
  }
  prediction <- offset + as.numeric(fixed %*% object$b) + as.numeric(mean)
  if(includeRandom){
    if(is.null(curveRandom)) prediction <- prediction+.nonlinear_random_prediction(object, newdata)
    else if(curveRandom$start) prediction <- prediction+as.numeric(
      design[, seq_len(curveRandom$start), drop=FALSE] %*% object$additiveRandomMode)
  }
  prediction
}

vcov.mmes.nonlinear <- function(object, ..., type=c("mean","all","working","random","conditional","group")){
  .nonlinear_check_dots(...)
  type <- match.arg(type)
  if(type=="group" && is.null(object$nonlinearRandom))
    stop("Group-parameter covariance requires random curve parameters.",call.=FALSE)
  uncertainty <- .nonlinear_uncertainty(object,if(type %in% c("random","conditional","group")) 2L else 0L)
  if(type=="group" && is.null(uncertainty$groupParameterCovariance))
    uncertainty <- post_nonlinear(object,mode=2L)$nonlinearInference
  switch(type, mean=uncertainty$parameterCovariance[names(stats::coef(object)),
      names(stats::coef(object)),drop=FALSE],
    all=uncertainty$parameterCovariance,working=uncertainty$workingCovariance,
    random=uncertainty$randomPEV,conditional=uncertainty$conditionalModeCovariance,
    group=uncertainty$groupParameterCovariance)
}

confint.mmes.nonlinear <- function(object, parm, level=0.95, ...,
    type=c("mean","all","working")){
  .nonlinear_check_dots(...)
  type <- match.arg(type)
  if(length(level)!=1L || !is.finite(level) || level<=0 || level>=1)
    stop("level must be strictly between zero and one.",call.=FALSE)
  estimate <- stats::coef(object,type=if(type=="mean") "mean" else type)
  if(missing(parm)) parm <- names(estimate)
  if(is.numeric(parm)){
    if(any(!is.finite(parm)) || any(parm!=floor(parm)) ||
       any(parm<1 | parm>length(estimate)))
      stop("Numeric parm must contain valid integer coefficient indices.",call.=FALSE)
    parm <- names(estimate)[parm]
  }
  if(!is.character(parm) || !length(parm) || anyNA(parm) ||
     any(!parm %in% names(estimate)) || anyDuplicated(parm))
    stop("parm must select unique coefficient names for the requested type.",call.=FALSE)
  uncertainty <- .nonlinear_uncertainty(object)
  covariance <- if(type=="working") uncertainty$workingCovariance else
    uncertainty$parameterCovariance
  center <- estimate[parm]
  error <- stats::qnorm((1+level)/2)*sqrt(diag(covariance)[parm])
  result <- cbind(center-error,center+error)
  if(type!="working"){
    positive <- paste0("curve:",object$nonlinear$parameters[object$nonlinear$positive])
    rows <- which(parm %in% positive)
    working <- uncertainty$workingCovariance
    center <- object$optimizer$par[parm[rows]]
    error <- stats::qnorm((1+level)/2)*sqrt(diag(working)[parm[rows]])
    result[rows,] <- exp(cbind(center-error,center+error))
  }
  colnames(result) <- paste0(format(100*c((1-level)/2,(1+level)/2),trim=TRUE),"%")
  result
}

anova.mmes.nonlinear <- function(object, object2=NULL, ..., parm=NULL, null=NULL){
  .nonlinear_check_dots(...)
  if(!is.null(object2)) return(.nonlinear_compare(object,object2))
  .nonlinear_wald(object,parm,null)
}

AIC.mmes.nonlinear <- function(object, ..., k=2){
  .nonlinear_information_criterion(object,list(...),k)
}

BIC.mmes.nonlinear <- function(object, ...){
  .nonlinear_information_criterion(object,list(...),NULL)
}

simulate.mmes.nonlinear <- function(object, nsim=1, seed=NULL,
    re.form=c("conditional","population"), maxMemoryMB=NULL, ...){
  .nonlinear_check_dots(...)
  .nonlinear_simulate(object,nsim,seed,match.arg(re.form),maxMemoryMB)
}

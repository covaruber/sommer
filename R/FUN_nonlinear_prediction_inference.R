.nonlinear_inference_proxy <- function(object){
  object$randomBlockWidths <- .nonlinear_block_widths(object)
  if(isTRUE(object$nonlinearFamily)) return(object)
  object$singleCategorical <- FALSE
  object$familySpecifications <- .laplace_family_specs(stats::gaussian(),list())$specifications
  object$familySpecifications[[1L]]$logParameter <- log(object$dispersion)
  object$laplaceFamilies <- list(stats::gaussian())
  object$laplaceFamilyIndex <- rep(1L,length(object$input))
  object$priorWeights <- rep(1,length(object$input))
  object
}

.nonlinear_prediction_shape <- function(values, template){
  if(is.list(template)){
    starts <- cumsum(c(0L,head(lengths(template),-1L)))
    result <- lapply(seq_along(template),function(index)
      .nonlinear_prediction_shape(values[starts[index]+seq_len(length(template[[index]]))],
        template[[index]]))
    names(result) <- names(template)
    return(result)
  }
  result <- template
  result[] <- values
  result
}

.nonlinear_observation_variance <- function(mean,specification,family,trials=1,rows=NULL){
  if(specification$family=="ziPoisson"){
    probability <- .laplace_zero_probability(specification,rows)
    return(mean+probability/(1-probability)*mean^2)
  }
  value <- family$variance(mean)
  if(specification$family=="binomial") value <- value/trials
  if(specification$family %in% c("gaussian","Gamma"))
    value <- value*exp(specification$logParameter)
  value
}

.nonlinear_predict_uncertainty <- function(object,newdata,re.form,type,
    se.fit,cov.fit,interval,level,trials){
  .mme_flag(se.fit,"se.fit")
  .mme_flag(cov.fit,"cov.fit")
  if(length(level)!=1L || !is.finite(level) || level<=0 || level>=1)
    stop("level must be strictly between zero and one.",call.=FALSE)
  conditional <- is.null(re.form)
  population <- identical(re.form,"population")
  uncertainty <- .nonlinear_uncertainty(object,if(conditional) 1L else 0L)
  numerical <- uncertainty$settings$numerical
  memory <- uncertainty$settings$memory
  parameters <- object$optimizer$par
  proxy <- .nonlinear_inference_proxy(object)
  evaluate <- function(parameters,mode=object$randomMode){
    local <- .nonlinear_inference_proxy(.nonlinear_at(object,parameters,mode))
    .nonlinear_family_predict(local,newdata,re.form,type)
  }
  point <- evaluate(parameters)
  count <- length(unlist(point,use.names=FALSE))
  .nonlinear_workspace(12*as.double(count)*(length(parameters)+
    if(conditional) length(object$randomMode) else 0)+
    if(cov.fit) 4*as.double(count)^2 else 0,memory,
    "Nonlinear prediction uncertainty workspace")
  derivative <- .nonlinear_jacobian(function(at){
    mode <- if(conditional) object$randomMode+
      as.numeric(uncertainty$modeSensitivity%*%(at-parameters)) else object$randomMode
    unlist(evaluate(at,mode),use.names=FALSE)
  },parameters,numerical,"Prediction parameter Jacobian")
  projected <- derivative%*%uncertainty$workingCovariance
  variance <- rowSums(projected*derivative)
  covariance <- if(cov.fit) projected%*%t(derivative) else NULL
  if(conditional){
    random <- .nonlinear_jacobian(function(mode)
      unlist(evaluate(parameters,mode),use.names=FALSE),
      object$randomMode,numerical,"Prediction random-effect Jacobian")
    projected <- random%*%uncertainty$conditionalModeCovariance
    variance <- variance+rowSums(projected*random)
    if(cov.fit) covariance <- covariance+projected%*%t(random)
  }
  if(any(!is.finite(variance)) || any(variance< -1e-10))
    stop("Nonlinear prediction variance is non-finite or negative.",call.=FALSE)
  standardError <- sqrt(pmax(0,variance))
  result <- list(fit=point,se.fit=.nonlinear_prediction_shape(standardError,point),
    method=paste(uncertainty$method,
      if(conditional) "; includes conditional random-effect curvature" else
        "; marginal parameter uncertainty only"))
  if(cov.fit) result$cov.fit <- covariance
  if(interval!="none"){
    intervalVariance <- variance
    if(interval=="prediction"){
      if(type!="response" || !is.null(proxy$categorical))
        stop("Future-observation prediction intervals require scalar response-scale predictions; categorical probabilities support confidence intervals only.",
          call.=FALSE)
      design <- .nonlinear_family_design(proxy,newdata,conditional || population)
      proxy$familySpecifications <- .laplace_prediction_specs(proxy,design$data,design$index)
      proxy$laplaceFamilies <- .laplace_update_families(proxy$laplaceFamilies,
        proxy$familySpecifications,design$index)
      n <- length(design$input)
      if(is.null(trials)){
        if(is.null(newdata)) trials <- proxy$priorWeights else{
          grouped <- any(vapply(proxy$familySpecifications,function(s) s$family=="binomial",
            logical(1))) && any(proxy$priorWeights>1)
          if(grouped)
            stop("Grouped-binomial prediction intervals for newdata require trials.",call.=FALSE)
          trials <- rep(1,n)
        }
      }
      if(!is.numeric(trials) || !length(trials) || !length(trials) %in% c(1L,n) ||
         any(!is.finite(trials)) || any(trials<1) || any(trials!=floor(trials)))
        stop("trials must contain positive integers, one or one per predicted response.",call.=FALSE)
      trials <- rep(trials,length.out=n)
      mean <- as.numeric(point)
      observationVariance <- numeric(n)
      if(population){
        zero <- .nonlinear_family_eta(proxy,design)
        second <- .nonlinear_family_population(proxy,design,zero,"response",
          moment=2L,trials=trials)
        observationVariance <- as.numeric(second)-mean^2
      }else for(index in seq_along(proxy$familySpecifications)){
        rows <- which(design$index==index)
        if(length(rows))
          observationVariance[rows] <- .nonlinear_observation_variance(mean[rows],
            proxy$familySpecifications[[index]],proxy$laplaceFamilies[[index]],
            trials[rows],rows)
      }
      if(any(!is.finite(observationVariance)) || any(observationVariance< -1e-8))
        stop("Future-response variance is non-finite or negative.",call.=FALSE)
      intervalVariance <- variance+pmax(0,observationVariance)
      result$method <- paste(result$method,
        "; normal future-observation approximation (not an exact discrete interval)")
    }
    quantile <- stats::qnorm((1+level)/2)
    result$lower <- .nonlinear_prediction_shape(as.numeric(unlist(point,use.names=FALSE))-
      quantile*sqrt(pmax(0,intervalVariance)),point)
    result$upper <- .nonlinear_prediction_shape(as.numeric(unlist(point,use.names=FALSE))+
      quantile*sqrt(pmax(0,intervalVariance)),point)
    result$level <- level
    result$interval <- interval
  }
  result
}

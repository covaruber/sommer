.spec_split_mode <- function(mode, widths){
  starts <- cumsum(c(0L,head(widths,-1L)))
  lapply(seq_along(widths),function(i) mode[starts[i]+seq_len(widths[i])])
}

.spec_conditional <- function(parameters, context, covariances){
  precision <- chol2inv(chol(as.matrix(Matrix::bdiag(covariances))))
  specifications <- .spec_family_proposal(context,parameters)
  previous <- NULL
  cached <- NULL
  evaluate <- function(mode){
    if(!identical(mode,previous)){
      predictor <- .spec_predictor_derivatives(parameters,mode,context)
      family <- .spec_family_evaluation(context,predictor$eta,specifications)
      cached <<- list(predictor=predictor,family=family)
      previous <<- mode
    }
    cached
  }
  value <- function(mode){
    evaluation <- evaluate(mode)
    -evaluation$family$logDensity+
      .5*sum(mode*(precision %*% mode))
  }
  gradient <- function(mode){
    evaluation <- evaluate(mode)
    as.numeric(-crossprod(evaluation$predictor$jacobian,evaluation$family$score)+precision %*% mode)
  }
  hessian <- function(mode){
    evaluation <- evaluate(mode)
    .spec_conditional_hessian(evaluation$predictor,evaluation$family,precision,context)
  }
  list(value=value,gradient=gradient,hessian=hessian)
}

.spec_future_variance <- function(object,newdata,target,trials){
  if(!is.null(object$baseFit)){
    mean <- predict(object,newdata,type="response",target=target)
    if(object$family$family=="binomial"){
      if(is.null(trials) && !is.null(newdata) && any(object$baseFit$prior.weights>1))
        stop("Grouped-binomial prediction intervals require trials for newdata.",call.=FALSE)
      if(is.null(trials)) trials <- if(is.null(newdata)) object$baseFit$prior.weights else 1
      if(!is.numeric(trials) || !length(trials) %in% c(1L,length(mean)) ||
         any(!is.finite(trials)) || any(trials<1 | trials!=floor(trials)))
        stop("trials must be positive integers, one or one per prediction row.",call.=FALSE)
      return(mean*(1-mean)/rep_len(trials,length(mean)))
    }
    return(object$family$variance(mean)*object$dispersion)
  }
  data <- .spec_prediction_data(object,newdata)
  mean <- predict(object,newdata,type="response",target=target)
  index <- .laplace_family_index(object$predictionFamily,data,nrow(data))
  binomial <- vapply(object$familySpecifications,function(s) s$family=="binomial",logical(1))
  if(is.null(trials)){
    grouped <- any(object$likelihoodContext$trials>1)
    if(!is.null(newdata) && grouped && any(binomial[index]))
      stop("Grouped-binomial prediction intervals require trials for newdata.",call.=FALSE)
    trials <- if(is.null(newdata)) object$likelihoodContext$trials else rep(1,nrow(data))
  }
  if(!is.numeric(trials) || !length(trials) %in% c(1L,nrow(data)) ||
     any(!is.finite(trials)) || any(trials<1 | trials!=floor(trials)))
    stop("trials must be positive integers, one or one per prediction row.",call.=FALSE)
  trials <- rep_len(trials,nrow(data))
  if(any(trials[!binomial[index]]!=1))
    stop("trials may be supplied only for binomial responses.",call.=FALSE)
  proxy <- list(familySpecifications=object$familySpecifications)
  specifications <- .laplace_prediction_specs(proxy,data,index)
  if(target=="population" && length(object$randomRecords)){
    fixed <- .spec_matrix(NULL,data,object$fixedTemplate)
    fixed$x <- .spec_fixed_prediction(object,data)
    design <- list(fixed=fixed,parameters=lapply(object$design$parameters,
      function(template) .spec_matrix(NULL,data,template)))
    random <- lapply(object$randomRecords,.spec_random_design,newdata=data,
      training=object$trainingData,control=.mmes_fitting_control(object))
    second <- .spec_population(object,data,design,random,"response",specifications,
      moment=2L,trials=trials)
    return(pmax(0,second-mean^2))
  }
  families <- .laplace_update_families(object$laplaceFamilies,specifications,index)
  result <- numeric(length(mean))
  for(position in seq_along(specifications)){
    rows <- which(index==position)
    result[rows] <- .nonlinear_observation_variance(mean[rows],specifications[[position]],
      families[[position]],trials[rows],rows)
  }
  result
}

.spec_report_parameters <- function(object, parameters){
  context <- object$likelihoodContext
  result <- stats::setNames(parameters[seq_len(context$p)],names(object$b))
  for(record in seq_along(context$records)){
    mapping <- context$records[[record]]$mapping
    working <- as.numeric(mapping$map$offset+mapping$map$T %*%
      parameters[context$covariancePositions[[record]]])
    start <- 0L
    for(term in seq_along(mapping$structures)){
      descriptor <- mapping$structures[[term]]
      values <- working[start+seq_along(descriptor$par)]
      start <- start+length(values)
      sigma <- as.matrix(evaluate_covstruct_cpp(descriptor,values))
      selected <- which(upper.tri(sigma,diag=TRUE),arr.ind=TRUE)
      result <- c(result,stats::setNames(sigma[selected],paste0("covariance",record,
        ":",term,":",selected[,1L],".",selected[,2L])))
    }
  }
  specifications <- .spec_family_proposal(context,parameters)
  for(position in seq_along(specifications)){
    specification <- specifications[[position]]
    if(is.null(specification$parameter)) next
    values <- if(specification$family=="ordinal") specification$thresholds else
      if(!is.null(specification$zeroDesign)) specification$zeroCoefficients else
        if(specification$family=="ziPoisson") stats::plogis(specification$logParameter) else
          exp(specification$logParameter)
    result <- c(result,stats::setNames(values,paste0("family:",
      names(specifications)[position],":",specification$parameter,".",seq_along(values))))
  }
  result
}

.spec_group_parameters <- function(object, parameters, mode){
  context <- object$likelihoodContext
  mode <- .spec_split_mode(mode,context$widths)
  values <- lapply(names(context$spec$parameters),function(name){
    declaration <- context$spec$parameters[[name]]
    value <- as.numeric(context$design$parameters[[name]]$x %*%
      parameters[declaration$positions])+context$design$parameters[[name]]$offset
    if(length(declaration$randomIndex)){
      record <- declaration$randomIndex
      value <- value+as.numeric(context$randomDesign[[record]] %*% mode[[record]])
    }
    if(declaration$link=="log") exp(value) else value
  })
  if(!length(values)) return(matrix(numeric(),context$n,0L))
  result <- do.call(cbind,values)
  colnames(result) <- names(context$spec$parameters)
  result
}

.spec_post_nonlinear <- function(object, mode=1L, control=NULL){
  if(!is.numeric(mode) || length(mode)!=1L || is.na(mode) || !mode %in% 0:2)
    stop("mode must be 0, 1 or 2.",call.=FALSE)
  if(!isTRUE(object$convergence))
    stop("Specification inference requires a converged fit.",call.=FALSE)
  if(!is.null(object$baseFit)){
    object$nonlinearInference <- list(meanCovariance=stats::vcov(object),
      workingCovariance=stats::vcov(object),mode=mode,
      method="Exact fixed-only model covariance")
    return(object)
  }
  settings <- .nonlinear_inference_settings(object,control)
  numerical <- settings$numerical
  context <- object$likelihoodContext
  context$control$task$nonlinear$mode <- .nonlinear_inference_mode_settings(
    context$control$task$nonlinear$mode)
  parameters <- object$optimizer$par
  q <- context$q
  .nonlinear_workspace(8*as.double(length(parameters)+q+context$n)^2,
    settings$memory,"Specification inference workspace")
  evaluate <- function(at) .spec_evaluate(at,context,TRUE)
  objective <- function(at) .spec_evaluate(at,context)
  accepted <- evaluate(parameters)
  if(abs(accepted$value+object$llik)>1e-7*(1+abs(object$llik)))
    stop("Reconstructed specification likelihood differs from the saved fit.",call.=FALSE)
  gradient <- .spec_gradient(parameters,objective)
  if(any(!is.finite(gradient)))
    stop("Specification outer score is non-finite.",call.=FALSE)
  if(accepted$competingModes>0L)
    stop("Competing conditional modes make local inference unavailable.",call.=FALSE)
  covarianceCoordinates <- unlist(context$covariancePositions,use.names=FALSE)
  nuisanceCoordinates <- seq.int(context$p+context$covarianceCount+1L,length(parameters))
  if(context$p+context$covarianceCount==length(parameters)) nuisanceCoordinates <- integer()
  bounded <- c(covarianceCoordinates,nuisanceCoordinates)
  if(length(bounded) && any(abs(parameters[bounded])>=25))
    stop("Specification covariance or nuisance coordinates are on a numerical boundary.",call.=FALSE)
  for(covariance in accepted$covariances){
    diagonal <- sqrt(diag(covariance))
    if(any(diagonal<exp(-12.5)) ||
       min(eigen(covariance/outer(diagonal,diagonal),symmetric=TRUE,
         only.values=TRUE)$values)<=numerical$curvatureTol)
      stop("Specification covariance is on a boundary.",call.=FALSE)
  }
  information <- .nonlinear_jacobian(function(at) .spec_gradient(at,objective),
    parameters,numerical,"Specification observed information")
  if(length(parameters)){
    working <- .nonlinear_positive_inverse(information,numerical,
      "Specification observed information")
    informationCondition <- attr(working,"condition")
    attr(working,"condition") <- NULL
    if(max(abs(gradient)*sqrt(diag(working)))>.1)
      stop("Specification outer score is not stationary on the estimated uncertainty scale.",call.=FALSE)
  }else{
    working <- information
    informationCondition <- 1
  }
  dimnames(working) <- list(names(parameters),names(parameters))
  reporting <- function(at) .spec_report_parameters(object,at)
  report <- reporting(parameters)
  derivative <- .nonlinear_jacobian(reporting,parameters,numerical,
    "Specification natural-parameter derivative")
  natural <- derivative %*% working %*% t(derivative)
  dimnames(natural) <- list(names(report),names(report))
  count <- seq_len(context$p)
  inference <- list(mode=as.integer(mode),workingCovariance=working,
    meanCovariance=working[count,count,drop=FALSE],parameterCovariance=natural,
    parameters=report,information=information,informationCondition=informationCondition,
    score=gradient,modeSettings=context$control$task$nonlinear$mode,
    method=paste(object$likelihoodMethod,"checked observed-information delta approximation"),
    randomMethod="Local empirical Bayes approximation; not an exact frequentist PEV")
  if(q && mode>0L){
    atMode <- unlist(accepted$mode,use.names=FALSE)
    conditional <- .spec_conditional(parameters,context,accepted$covariances)
    numericalHessian <- .nonlinear_jacobian(conditional$gradient,atMode,numerical,
      "Specification conditional Hessian")
    hessian <- conditional$hessian(atMode)
    if(any(abs(hessian-numericalHessian)>
      numerical$derivativeTol*(1+abs(hessian))))
      stop("Specification conditional curvature disagrees with mode-score derivatives.",
        call.=FALSE)
    conditionalCovariance <- .nonlinear_positive_inverse(hessian,numerical,
      "Specification conditional Hessian")
    attr(conditionalCovariance,"condition") <- NULL
    sensitivity <- .nonlinear_jacobian(function(at)
      unlist(evaluate(at)$mode,use.names=FALSE),parameters,numerical,
      "Specification conditional-mode sensitivity")
    pev <- conditionalCovariance+sensitivity %*% working %*% t(sensitivity)
    inference$modeSensitivity <- sensitivity
    inference$conditionalCovariance <- conditionalCovariance
    inference$randomPEV <- if(mode==2L) pev else diag(pev)
    if(mode==2L) inference$parameterModeCovariance <- working %*% t(sensitivity)
    group <- .spec_group_parameters(object,parameters,atMode)
    if(ncol(group)){
      .nonlinear_workspace(as.double(length(group))*
        (4*length(parameters)+4*q+if(mode==2L) 4*length(group) else 4),
        settings$memory,"Specification group-parameter workspace")
      total <- .nonlinear_jacobian(function(at)
        .spec_group_parameters(object,at,unlist(evaluate(at)$mode,use.names=FALSE)),
        parameters,numerical,"Specification group-parameter derivative")
      random <- .nonlinear_jacobian(function(b)
        .spec_group_parameters(object,parameters,b),atMode,numerical,
        "Specification group-mode derivative")
      variance <- rowSums((total %*% working)*total)+
        rowSums((random %*% conditionalCovariance)*random)
      inference$groupParameterSE <- matrix(sqrt(pmax(0,variance)),
        nrow(group),ncol(group),dimnames=dimnames(group))
      if(mode==2L) inference$groupCovariance <- total %*% working %*% t(total)+
        random %*% conditionalCovariance %*% t(random)
    }
  }
  object$nonlinearInference <- inference
  if(q && mode>0L){
    object$uPevList <- .spec_split_mode(if(mode==2L) diag(inference$randomPEV) else
      inference$randomPEV,context$widths)
    object$uPevMode <- as.integer(mode)
  }
  object
}

vcDelta.mmes.specification <- function(object, transform){
  if(!inherits(transform,"formula") || !length(transform) %in% 2:3)
    stop("transform must be a one- or two-sided formula.",call.=FALSE)
  if(!is.null(object$baseFit))
    stop("Variance-component delta transforms require a specification covariance layout.",
      call.=FALSE)
  post <- .spec_post_nonlinear(object,0L)
  estimate <- post$nonlinearInference$parameters
  selected <- seq_along(estimate)>length(object$b)
  estimate <- estimate[selected]
  if(!length(estimate)) stop("No covariance or family parameters to transform.",call.=FALSE)
  covariance <- post$nonlinearInference$parameterCovariance[selected,selected,drop=FALSE]
  parameters <- stats::setNames(as.list(estimate),paste0("V",seq_along(estimate)))
  expression <- transform[[length(transform)]]
  value <- eval(deriv(expression,names(parameters)),parameters)
  gradient <- as.numeric(attr(value,"gradient"))
  if(length(value)!=1L || any(!is.finite(c(value,gradient))))
    stop("Variance transform and its derivative must be finite scalar values.",call.=FALSE)
  variance <- as.numeric(crossprod(gradient,covariance %*% gradient))
  if(variance< -1e-10)
    stop("Variance transform has negative propagated variance.",call.=FALSE)
  result <- data.frame(Estimate=as.numeric(value),SE=sqrt(max(0,variance)))
  attr(result,"parameters") <- stats::setNames(names(estimate),names(parameters))
  class(result) <- c("vcDelta.mmes","data.frame")
  result
}

confint.mmes.specification <- function(object, parm, level=.95, ...,
    type=c("mean","all","working")){
  .spec_dots(...)
  type <- match.arg(type)
  if(!is.numeric(level) || length(level)!=1L || !is.finite(level) || level<=0 || level>=1)
    stop("level must be strictly between zero and one.",call.=FALSE)
  if(!is.null(object$baseFit)) return(stats::confint.default(object$baseFit,
    parm=if(missing(parm)) seq_along(object$b) else parm,level=level))
  post <- .spec_post_nonlinear(object,0L)
  estimate <- if(type=="all") post$nonlinearInference$parameters else
    if(type=="working") object$optimizer$par else object$b
  covariance <- stats::vcov(post,type=type)
  if(missing(parm)) parm <- names(estimate)
  if(is.numeric(parm)) parm <- names(estimate)[parm]
  if(anyNA(parm) || anyDuplicated(parm) || any(!parm %in% names(estimate)))
    stop("parm must select unique known parameters.",call.=FALSE)
  margin <- stats::qnorm((1+level)/2)*sqrt(diag(covariance)[parm])
  result <- cbind(estimate[parm]-margin,estimate[parm]+margin)
  colnames(result) <- paste0(format(100*c((1-level)/2,(1+level)/2)),"%")
  result
}

.spec_prediction_inference <- function(object,newdata,type,target,se.fit,interval,level,
    cov.fit=FALSE,trials=NULL){
  if(!is.numeric(level) || length(level)!=1L || !is.finite(level) || level<=0 || level>=1)
    stop("level must be strictly between zero and one.",call.=FALSE)
  if(!is.null(object$baseFit)){
    result <- if(inherits(object$baseFit,"glm"))
      stats::predict(object$baseFit,newdata=newdata,type=type,se.fit=TRUE,
        dispersion=object$dispersion) else
      stats::predict(object$baseFit,newdata=newdata,se.fit=TRUE,
        scale=sqrt(object$dispersion))
    fit <- result$fit
    design <- if(is.null(newdata)) stats::model.matrix(object$baseFit) else
      stats::model.matrix(stats::delete.response(stats::terms(object$baseFit)),
        data=newdata,contrasts.arg=object$baseFit$contrasts,xlev=object$baseFit$xlevels)
    if(nrow(design)!=length(fit) || any(!is.finite(design)))
      stop("Prediction design has missing or nonfinite rows.",call.=FALSE)
    if(type=="response" && inherits(object$baseFit,"glm")){
      eta <- if(is.null(newdata)) object$baseFit$linear.predictors else
        stats::predict(object$baseFit,newdata=newdata,type="link")
      design <- design*as.numeric(object$family$mu.eta(eta))
    }
    covariance <- design %*% stats::vcov(object) %*% t(design)
  }else{
    post <- .spec_post_nonlinear(object,if(target=="conditional") 1L else 0L)
    inference <- post$nonlinearInference
    parameters <- object$optimizer$par
    settings <- .nonlinear_inference_settings(object)$numerical
    fit <- predict(object,newdata,type=type,target=target)
    restore <- function(at){
      result <- object
      evaluation <- .spec_evaluate(at,object$likelihoodContext,TRUE)
      result$b <- at[seq_along(object$b)]
      result$optimizer$par <- at
      result$randomMode <- evaluation$mode
      result$randomCovariances <- evaluation$covariances
      result$familySpecifications <- evaluation$specifications
      result$fitted.values <- .spec_response(result,evaluation$eta,"response",result$trainingData)
      result$linear.predictors <- evaluation$eta
      result
    }
    derivative <- .nonlinear_jacobian(function(at)
      predict(restore(at),newdata,type=type,target=target),parameters,settings,
      "Specification prediction derivative")
    covariance <- derivative %*% inference$workingCovariance %*% t(derivative)
    if(target=="conditional" && object$likelihoodContext$q){
      random <- unlist(object$randomMode,use.names=FALSE)
      derivative <- .nonlinear_jacobian(function(mode){
        result <- object
        result$randomMode <- .spec_split_mode(mode,object$likelihoodContext$widths)
        predict(result,newdata,type=type,target=target)
      },random,settings,"Specification prediction-mode derivative")
      covariance <- covariance+
        derivative %*% inference$conditionalCovariance %*% t(derivative)
    }
  }
  uncertainty <- .nonlinear_prediction_shape(sqrt(pmax(0,diag(covariance))),fit)
  result <- list(fit=fit,se.fit=uncertainty,
    approximation="Local delta approximation; conditional predictions include empirical Bayes mode uncertainty")
  if(cov.fit) result$cov.fit <- covariance
  if(interval!="none"){
    variance <- diag(covariance)
    if(interval=="prediction"){
      if(type!="response" || is.list(fit) || is.matrix(fit))
        stop("Prediction intervals require scalar response-scale outputs.",call.=FALSE)
      variance <- variance+.spec_future_variance(object,newdata,target,trials)
    }
    margin <- stats::qnorm((1+level)/2)*sqrt(pmax(0,variance))
    result$lower <- .nonlinear_prediction_shape(as.numeric(unlist(fit))-margin,fit)
    result$upper <- .nonlinear_prediction_shape(as.numeric(unlist(fit))+margin,fit)
  }
  result
}

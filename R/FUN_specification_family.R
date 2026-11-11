.spec_family_prepare <- function(fixed, data, family, familyParameters, control,
    allowLegacyLinks=FALSE){
  reserved <- ".specificationRow"
  if(reserved %in% names(data))
    stop("Specification fitting reserves .specificationRow.",call.=FALSE)
  data[[reserved]] <- seq_len(nrow(data))
  local <- control
  local$task$glmm$familyParameters <- familyParameters
  preparationFamily <- family
  adapter <- NULL
  if(allowLegacyLinks && inherits(family,"family") &&
     family$family %in% c("gaussian","poisson","binomial")){
    links <- switch(family$family,gaussian="identity",poisson="log",
      binomial=c("logit","probit","cloglog"))
    if(!family$link %in% links){
      preparationFamily <- switch(family$family,gaussian=stats::gaussian(),
        poisson=stats::poisson(),binomial=stats::binomial())
      preparationFamily$trials <- family$trials
      adapter <- list(family=family,canonical=preparationFamily)
    }
  }
  prepared <- .nonlinear_family_prepare(fixed,NULL,data,preparationFamily,reserved,local,
    "exclude","exclude","exclude","exclude",NULL)
  retained <- prepared$data[prepared$included,,drop=FALSE]
  metadata <- prepared$prepared$categorical
  prepared$categoricalModel <- !is.null(metadata)
  prepared$predictionFamily <- prepared$family
  prepared$linkAdapter <- adapter
  if(!is.null(adapter)){
    prepared$families[[1L]] <- family
    prepared$predictionFamily <- family
  }
  if(!is.null(metadata)){
    prepared$fixedMap <- metadata$map
    prepared$fixedBasis <- metadata$basisColumns
  }
  prepared$specifications <- .laplace_zero_designs(prepared$families,
    prepared$specifications,retained,rep(TRUE,nrow(retained)),prepared$index)
  prepared
}

.spec_family_proposal <- function(context, at){
  .laplace_proposal_specs(context$familySpecifications,at,
    list(fixedCount=context$p,covarianceCount=context$covarianceCount))
}

.spec_family_evaluation <- function(context, eta, specifications){
  adapter <- context$linkAdapter
  if(!is.null(adapter)){
    mean <- adapter$family$linkinv(eta)
    predictor <- adapter$canonical$linkfun(mean)
    jacobian <- adapter$family$mu.eta(eta)/adapter$canonical$mu.eta(predictor)
  }else predictor <- eta
  result <- sommer_glmm_family_eval(context$y,as.numeric(predictor),specifications,
    context$familyIndex,context$trials,memoryMB=context$control$execution$denseMemoryMB)
  if(!is.null(adapter)){
    derivative <- function(at){
      predictor <- adapter$canonical$linkfun(adapter$family$linkinv(at))
      adapter$family$mu.eta(at)/adapter$canonical$mu.eta(predictor)
    }
    step <- .Machine$double.eps^(1/3)*(1+abs(eta))
    second <- (derivative(eta+step)-derivative(eta-step))/(2*step)
    result$curvature <- as.numeric(result$curvature)*jacobian^2-
      as.numeric(result$score)*second
    result$score <- as.numeric(result$score)*jacobian
    result$blockCurvature <- NULL
  }
  result
}

.spec_family_proxy <- function(object, data, specifications=object$familySpecifications){
  index <- .laplace_family_index(object$predictionFamily,data,nrow(data))
  rowMap <- attr(data,"specificationRowMap")
  coordinate <- attr(data,"specificationCoordinate")
  if(is.null(rowMap)) rowMap <- seq_len(nrow(data))
  if(is.null(coordinate)) coordinate <- rep(1L,nrow(data))
  list(familySpecifications=specifications,laplaceFamilyIndex=index,
    laplaceFamilies=.laplace_update_families(object$laplaceFamilies,specifications,index),
    laplacePostPEVSettings=list(denseMemoryMB=object$control$execution$denseMemoryMB),
    categorical=list(rowMap=rowMap,coordinate=coordinate),singleCategorical=FALSE)
}

.spec_response <- function(object, eta, type, data, specifications=object$familySpecifications){
  if(any(!is.finite(eta)))
    stop("Specification prediction produced nonfinite link predictors.",call.=FALSE)
  if(!object$categoricalModel){
    if(type=="link") result <- as.numeric(eta) else{
      index <- .laplace_family_index(object$predictionFamily,data,nrow(data))
      families <- .laplace_update_families(object$laplaceFamilies,specifications,index)
      result <- .pql_apply(families,index,"linkinv",eta)
    }
  }else{
    proxy <- .spec_family_proxy(object,data,specifications)
    result <- .laplace_mixed_output(proxy,eta,type)
    if(object$singleCategorical) result <- result[[1L]]
  }
  if(any(!is.finite(unlist(result,use.names=FALSE))))
    stop("Specification prediction produced nonfinite response values.",call.=FALSE)
  result
}

.spec_prediction_data <- function(object, newdata){
  if(is.null(newdata)) return(object$trainingData)
  if(!is.data.frame(newdata) || !nrow(newdata))
    stop("Specification prediction requires nonempty newdata.",call.=FALSE)
  width <- max(vapply(object$familySpecifications,function(s)
    if(s$family=="multinomial") length(s$predictorLevels) else 1L,integer(1)))
  modes <- sum(vapply(object$randomMode,length,integer(1)))
  .nonlinear_workspace(as.double(nrow(newdata))*width*
    (40+length(object$b)+modes),object$control$execution$denseMemoryMB,
    "Specification prediction workspace")
  if(!object$categoricalModel) return(newdata)
  if(object$singleCategorical) newdata$.nonlinearTrait <- "response"
  expanded <- .laplace_mixed_expand(newdata,object$predictionFamily,
    object$familySpecifications)
  attr(expanded$data,"specificationRowMap") <- expanded$rowMap
  attr(expanded$data,"specificationCoordinate") <- expanded$coordinate
  expanded$data
}

.spec_fixed_prediction <- function(object, data){
  basis <- .spec_matrix(NULL,data,object$fixedTemplate)$x
  if(is.null(object$fixedMap)) return(basis)
  colnames(basis)[colnames(basis)=="(Intercept)"] <- "Intercept"
  basis <- basis[,object$fixedBasis,drop=FALSE]
  proxy <- .spec_family_proxy(object,data)
  as.matrix(.laplace_mixed_design(basis,proxy$laplaceFamilyIndex,
    proxy$categorical$coordinate,object$fixedMap))
}

.spec_population <- function(object, data, design, random, type, specifications,
    moment=1L,trials=NULL){
  index <- .laplace_family_index(object$predictionFamily,data,nrow(data))
  mode <- lapply(object$randomMode,function(b) numeric(length(b)))
  zero <- .spec_predictor(object$specification,object$b,data,design,random,mode)
  output <- .spec_response(object,zero,type,data,specifications)
  if(!object$categoricalModel || object$singleCategorical) output <- list(output)
  settings <- .mmes_fitting_control(object)$task$nonlinear$population
  for(position in seq_along(specifications)){
    component <- specifications[[position]]
    rows <- which(index==position)
    if(!length(rows)) next
    d <- if(component$family=="multinomial") length(component$predictorLevels) else 1L
    groups <- matrix(rows,ncol=d,byrow=TRUE)
    values <- lapply(seq_len(nrow(groups)),function(i){
      current <- groups[i,]
      covariances <- lapply(seq_along(random),function(j){
        z <- random[[j]][current,,drop=FALSE]
        z %*% object$randomCovariances[[j]] %*% t(z)
      })
      covariance <- as.matrix(Matrix::bdiag(covariances))
      evaluate <- function(deviations){
        count <- nrow(deviations)
        selected <- rep(current,times=count)
        localDesign <- list(fixed=list(x=design$fixed$x[selected,,drop=FALSE],
          offset=design$fixed$offset[selected]),
          parameters=lapply(design$parameters,function(x)
            list(x=x$x[selected,,drop=FALSE],offset=x$offset[selected])))
        shifts <- lapply(seq_along(random),function(j)
          as.numeric(t(deviations[,(j-1L)*d+seq_len(d),drop=FALSE])))
        eta <- .spec_predictor(object$specification,object$b,data[selected,,drop=FALSE],
          localDesign,lapply(random,function(z) z[selected,,drop=FALSE]),mode,shifts)
        eta <- matrix(eta,count,d,byrow=TRUE)
        if(type=="link") return(eta)
        if(component$family=="multinomial")
          return(.laplace_nominal_probabilities(eta,component))
        if(component$family=="ordinal")
          return(.laplace_ordinal_distribution(as.numeric(eta),component,
            component$levels,object$control$execution$denseMemoryMB)$probability)
        family <- .laplace_update_families(object$laplaceFamilies,specifications)[[position]]
        mu <- if(component$family=="ziPoisson")
          (1-.laplace_zero_probability(component,current))*exp(eta) else family$linkinv(eta)
        if(moment==1L) return(mu)
        mu^2+.nonlinear_observation_variance(as.numeric(mu),component,family,
          if(is.null(trials)) 1 else trials[current],current)
      }
      .nonlinear_normal_integral(covariance,evaluate,settings,
        object$control$execution$denseMemoryMB)
    })
    value <- do.call(rbind,values)
    if(!object$categoricalModel) output[[1L]][rows] <- as.numeric(value) else
      output[[position]] <- if(is.matrix(output[[position]])){
        existing <- output[[position]]
        existing[] <- value
        existing
      }else as.numeric(value)
  }
  if(!object$categoricalModel || object$singleCategorical) output[[1L]] else output
}

.spec_draw_response <- function(object, eta, specifications){
  data <- object$trainingData
  index <- object$likelihoodContext$familyIndex
  trials <- object$likelihoodContext$trials
  responses <- object$likelihoodContext$y
  proxy <- .spec_family_proxy(object,data,specifications)
  for(position in seq_along(specifications)){
    rows <- which(index==position)
    specification <- specifications[[position]]
    family <- specification$family
    predictor <- eta[rows]
    if(family=="multinomial"){
      groups <- matrix(rows,ncol=length(specification$predictorLevels),byrow=TRUE)
      probabilities <- .laplace_nominal_probabilities(
        matrix(predictor,nrow=nrow(groups),byrow=TRUE),specification)
      draw <- vapply(seq_len(nrow(groups)),function(i)
        sample.int(ncol(probabilities),1L,prob=probabilities[i,]),integer(1))
      responses[rows] <- rep(draw,each=ncol(groups))
    }else if(family=="ordinal"){
      probability <- .laplace_ordinal_distribution(predictor,specification,
        specification$levels,object$control$execution$denseMemoryMB)$probability
      responses[rows] <- vapply(seq_along(rows),function(i)
        sample.int(ncol(probability),1L,prob=probability[i,]),integer(1))
    }else{
      mu <- proxy$laplaceFamilies[[position]]$linkinv(predictor)
      responses[rows] <- switch(family,
      gaussian=stats::rnorm(length(rows),mu,exp(specification$logParameter/2)),
      binomial=stats::rbinom(length(rows),trials[rows],mu)/trials[rows],
      poisson=stats::rpois(length(rows),mu),
      nbinom2=stats::rnbinom(length(rows),size=exp(specification$logParameter),mu=mu),
      ziPoisson={
        count <- stats::rpois(length(rows),exp(predictor))
        count[stats::runif(length(rows))<.laplace_zero_probability(specification,rows)] <- 0
        count
      },
      Gamma=stats::rgamma(length(rows),shape=exp(-specification$logParameter),
        scale=mu*exp(specification$logParameter)))
    }

  }
  if(any(!is.finite(responses)))
    stop("Specification simulation produced nonfinite responses.",call.=FALSE)
  if(!object$categoricalModel) return(responses)
  result <- lapply(seq_along(specifications),function(position){
    rows <- which(index==position)
    component <- specifications[[position]]
    if(component$family=="multinomial")
      rows <- rows[seq.int(1L,length(rows),by=length(component$predictorLevels))]
    value <- responses[rows]
    if(component$family %in% c("ordinal","multinomial"))
      factor(value,levels=seq_along(component$levels),labels=component$levels,
        ordered=component$family=="ordinal") else value
  })
  names(result) <- names(specifications)
  if(object$singleCategorical) result[[1L]] else result
}

.spec_categorical_residuals <- function(object, type){
  result <- lapply(seq_along(object$familySpecifications),function(position){
    component <- object$familySpecifications[[position]]
    rows <- which(object$likelihoodContext$familyIndex==position)
    observed <- object$y[rows]
    fitted <- if(object$singleCategorical) object$fitted.values else object$fitted.values[[position]]
    if(component$family=="multinomial"){
      d <- length(component$predictorLevels)
      observed <- observed[seq.int(1L,length(rows),by=d)]
      indicator <- outer(observed,seq_along(component$levels),"==")*1
      residual <- indicator-fitted
      if(type=="pearson") residual <- residual/sqrt(fitted*(1-fitted))
      return(residual)
    }
    if(component$family=="ordinal"){
      values <- seq_along(component$levels)
      mean <- as.numeric(fitted %*% values)
      variance <- as.numeric(fitted %*% values^2)-mean^2
    }else{
      mean <- fitted
      family <- .laplace_update_families(object$laplaceFamilies,
        object$familySpecifications)[[position]]
      variance <- .nonlinear_observation_variance(mean,component,family,
        object$likelihoodContext$trials[rows],rows)
    }
    residual <- observed-mean
    if(type=="pearson") residual <- residual/sqrt(variance)
    residual
  })
  names(result) <- names(object$familySpecifications)
  if(object$singleCategorical) result[[1L]] else result
}

.nonlinear_family_output <- function(object, eta, type="response", prediction=NULL){
  if(!is.null(object$categorical)){
    value <- if(is.null(prediction)) .laplace_mixed_output(object,eta,type) else
      .laplace_mixed_output(object,eta,type,index=prediction$index,
        rowMap=prediction$rowMap,coordinate=prediction$coordinate)
    if(object$singleCategorical) value <- value[[1L]]
  }else{
    index <- if(is.null(prediction)) object$laplaceFamilyIndex else prediction$index
    value <- if(type=="link") eta else .pql_apply(object$laplaceFamilies,index,"linkinv",eta)
  }
  if(any(!is.finite(unlist(value,use.names=FALSE))))
    stop("Nonlinear prediction produced non-finite response values.",call.=FALSE)
  value
}

.nonlinear_family_design <- function(object,newdata,random){
  if(is.null(newdata)){
    n <- length(object$input)
    return(list(fixed=object$fixedDesign,random=object$randomDesign,
      offset=object$offset,input=object$input,index=object$laplaceFamilyIndex,
      rowMap=if(is.null(object$categorical)) seq_len(n) else object$categorical$rowMap,
      coordinate=if(is.null(object$categorical)) rep(1L,n) else
        object$categorical$coordinate,data=NULL))
  }
  if(!is.data.frame(newdata) || !nrow(newdata))
    stop("Nonlinear prediction requires nonempty newdata in a data frame.",call.=FALSE)
  if(isTRUE(object$singleCategorical)){
    newdata$.nonlinearTrait <- "response"
  }
  if(8*as.double(nrow(newdata))*(40+ncol(object$fixedDesign)) >
     object$control$execution$denseMemoryMB*1048576)
    stop("Nonlinear prediction workspace exceeds execution$denseMemoryMB.",call.=FALSE)
  if(!is.null(object$categorical)){
    local <- object
    local$family <- object$predictionFamily
    local$predictionSettings$mmesArgs$observationMask <- NULL
    design <- .laplace_mixed_prediction_design(local,newdata,
      includeRandom=random,population=FALSE)
    expanded <- .laplace_mixed_expand(newdata,object$predictionFamily,object$familySpecifications)$data
    design$input <- .nonlinear_input(object$nonlinear,expanded)
    design$data <- expanded
    return(design)
  }
  input <- .nonlinear_input(object$nonlinear,newdata)
  frame <- stats::model.frame(object$fixedTerms,data=newdata,
    xlev=object$fixedXlevels,na.action=stats::na.pass)
  fixed <- stats::model.matrix(object$fixedTerms,frame,contrasts.arg=object$contrasts)
  colnames(fixed)[colnames(fixed)=="(Intercept)"] <- "Intercept"
  if(nrow(fixed)!=nrow(newdata) || !identical(colnames(fixed),object$predictionColumns) ||
     any(!is.finite(fixed)))
    stop("Prediction fixed design has missing values or differs from the fitted model.",call.=FALSE)
  offset <- stats::model.offset(frame)
  if(is.null(offset)) offset <- numeric(nrow(newdata))
  if(any(!is.finite(offset))) stop("Prediction offsets must be finite.",call.=FALSE)
  index <- .laplace_family_index(object$family,newdata,nrow(newdata))
  if(anyNA(index)) stop("Prediction family indices must not be missing.",call.=FALSE)
  z <- if(random && ncol(object$randomDesign)>0L)
    .nonlinear_random_prediction(object,newdata,design=TRUE) else
      Matrix::Matrix(0,nrow=nrow(newdata),ncol=0L,sparse=TRUE)
  list(fixed=fixed,random=z,offset=offset,input=input,index=index,
    rowMap=seq_len(nrow(newdata)),coordinate=rep(1L,nrow(newdata)),data=newdata)
}

.nonlinear_family_eta <- function(object,prediction,mode=NULL){
  layout <- object$nonlinearRandom
  n <- length(prediction$input)
  curve <- if(is.null(layout) || is.null(mode)){
    sommer_nonlinear_curve_eval(object$nonlinear$curve,prediction$input,
      object$nonlinearCoordinates,gradient=FALSE)$mean
  }else if(identical(layout$version,2L)){
    .nonlinear_curve_value(object,prediction$input,prediction$random,mode)
  }else{
    deviations <- as.numeric(prediction$random[,layout$start+seq_len(layout$count),
      drop=FALSE] %*% mode[layout$start+seq_len(layout$count)])
    sommer_nonlinear_curve_random_eval(object$nonlinear$curve,prediction$input,
      object$nonlinearCoordinates,layout$coordinate,deviations)
  }
  eta <- prediction$offset+as.numeric(prediction$fixed%*%object$b)+as.numeric(curve)
  if(!is.null(mode)){
    count <- if(is.null(layout)) length(mode) else layout$start
    if(count) eta <- eta+as.numeric(prediction$random[,seq_len(count),drop=FALSE]%*%mode[seq_len(count)])
  }
  if(length(eta)!=n || any(!is.finite(eta)))
    stop("Nonlinear prediction produced non-finite predictors.",call.=FALSE)
  eta
}

.nonlinear_family_predict <- function(object,newdata,re.form,type){
  if(!is.null(re.form) && !identical(re.form,NA) && !identical(re.form,NA_integer_) &&
     !identical(re.form,"population"))
    stop("re.form must be NULL, NA or 'population' for nonlinear models.",call.=FALSE)
  population <- identical(re.form,"population")
  design <- .nonlinear_family_design(object,newdata,is.null(re.form) || population)
  object$familySpecifications <- .laplace_prediction_specs(object,design$data,design$index)
  object$laplaceFamilies <- .laplace_update_families(object$laplaceFamilies,
    object$familySpecifications,design$index)
  eta <- .nonlinear_family_eta(object,design,if(is.null(re.form)) object$randomMode else NULL)
  if(population) return(.nonlinear_family_population(object,design,eta,type))
  .nonlinear_family_output(object,eta,type,design)
}

.nonlinear_normal_integral <- function(covariance,evaluate,settings,memoryMB){
  decomposition <- eigen((covariance+t(covariance))/2,symmetric=TRUE)
  tolerance <- 1e-12*max(1,max(abs(decomposition$values)))
  if(any(decomposition$values < -tolerance))
    stop("Population predictor covariance is not positive semidefinite.",call.=FALSE)
  selected <- which(decomposition$values>tolerance)
  if(!length(selected)){
    value <- as.numeric(evaluate(matrix(0,1L,nrow(covariance))))
    if(any(!is.finite(value)))
      stop("Nonlinear population prediction produced non-finite response values.",call.=FALSE)
    return(value)
  }
  loading <- sweep(decomposition$vectors[,selected,drop=FALSE],2L,
    sqrt(decomposition$values[selected]),"*")
  integrate <- function(order){
    rank <- length(selected)
    count <- as.double(order)^rank
    if(!is.finite(count) || count>.Machine$integer.max ||
       8*count*(8*rank+12*nrow(covariance)+24)>memoryMB*1048576)
      stop("Nonlinear population quadrature workspace exceeds execution$denseMemoryMB; reduce quadrature order or predictor dimension.",
        call.=FALSE)
    rule <- .laplace_normal_quadrature(order)
    grid <- as.matrix(expand.grid(rep(list(seq_len(order)),rank),
      KEEP.OUT.ATTRS=FALSE))
    nodes <- matrix(rule$nodes[grid],nrow(grid),rank)
    weights <- apply(matrix(rule$weights[grid],nrow(grid),rank),1L,prod)
    values <- as.matrix(evaluate(nodes%*%t(loading)))
    if(nrow(values)!=count || any(!is.finite(values)))
      stop("Nonlinear population quadrature produced non-finite response values.",call.=FALSE)
    as.numeric(crossprod(weights,values))
  }
  order <- settings$order
  previous <- integrate(order)
  repeat{
    order <- min(2L*order,settings$maxOrder)
    current <- integrate(order)
    if(all(abs(current-previous)<=settings$tol*(1+abs(current)))) return(current)
    if(order>=settings$maxOrder)
      stop("Nonlinear population quadrature did not meet task$nonlinear$population$tol at maxOrder.",
        call.=FALSE)
    previous <- current
  }
}

.nonlinear_family_population <- function(object,design,zero,type,moment=1L,trials=NULL){
  n <- length(zero)
  if(ncol(design$random)==0L){
    if(moment==1L) return(.nonlinear_family_output(object,zero,type,design))
    mean <- as.numeric(.nonlinear_family_output(object,zero,"response",design))
    variance <- numeric(n)
    for(index in seq_along(object$familySpecifications)){
      rows <- which(design$index==index)
      if(length(rows))
        variance[rows] <- .nonlinear_observation_variance(mean[rows],
          object$familySpecifications[[index]],object$laplaceFamilies[[index]],
          if(is.null(trials)) 1 else trials[rows],rows)
    }
    return(mean^2+variance)
  }
  layout <- object$nonlinearRandom
  additive <- design$random
  coordinates <- if(is.null(layout)) integer() else layout$coordinate
  bases <- .nonlinear_curve_bases(object,design$random)
  if(!is.null(layout)){
    columns <- layout$start+seq_len(layout$count)
    additive[,columns] <- 0
  }
  combined <- do.call(rbind,c(list(additive),bases[coordinates]))
  output <- .nonlinear_family_output(object,zero,type,design)
  if(object$singleCategorical) output <- list(output)
  if(is.null(object$categorical)) output <- list(response=output)
  for(position in seq_along(object$familySpecifications)){
    specification <- object$familySpecifications[[position]]
    rows <- which(design$index==position)
    if(!length(rows)) next
    d <- if(specification$family=="multinomial") length(specification$predictorLevels) else 1L
    groups <- matrix(rows,ncol=d,byrow=TRUE)
    covariance <- .laplace_latent_covariance(combined,object$covStruct,
      object$randomRelationships,object$randomBlockWidths,
      groups=do.call(cbind,lapply(0:length(coordinates),function(j) groups+j*n)),
      memoryMB=object$control$execution$denseMemoryMB,
      allowSingular=.nonlinear_has_reduced_rank(object$covStruct))
    values <- vector("list",nrow(groups))
    for(i in seq_len(nrow(groups))){
      current <- groups[i,]
      evaluate <- function(deviations){
        count <- nrow(deviations)
        eta <- matrix(rep(zero[current],each=count),count,d)+deviations[,seq_len(d),drop=FALSE]
        if(!is.null(layout)){
          shifts <- matrix(0,count*d,length(object$nonlinearCoordinates))
          for(j in seq_along(coordinates))
            shifts[,coordinates[j]] <- as.numeric(
              deviations[,j*d+seq_len(d),drop=FALSE])
          curve <- sommer_nonlinear_curve_joint_eval(object$nonlinear$curve,
            rep(design$input[current],each=count),object$nonlinearCoordinates,shifts)
          base <- sommer_nonlinear_curve_eval(object$nonlinear$curve,
            design$input[current],object$nonlinearCoordinates,gradient=FALSE)$mean
          eta <- eta+matrix(as.numeric(curve)-rep(as.numeric(base),each=count),count,d)
        }
        if(type=="link") return(eta)
        if(specification$family=="multinomial")
          return(.laplace_nominal_probabilities(eta,specification))
        if(specification$family=="ordinal")
          return(.laplace_ordinal_distribution(as.numeric(eta),specification,
            specification$levels,object$control$execution$denseMemoryMB)$probability)
        mean <- if(specification$family=="ziPoisson")
          (1-.laplace_zero_probability(specification,current))*exp(eta) else
            object$laplaceFamilies[[position]]$linkinv(eta)
        if(moment==1L) return(mean)
        mean^2+.nonlinear_observation_variance(as.numeric(mean),specification,
          object$laplaceFamilies[[position]],
          if(is.null(trials)) 1 else trials[current],current)
      }
      localCovariance <- matrix(covariance[,,i],nrow=dim(covariance)[1L])
      values[[i]] <- .nonlinear_normal_integral(localCovariance,evaluate,
        .mmes_fitting_control(object)$task$nonlinear$population,object$control$execution$denseMemoryMB)
    }
    value <- do.call(rbind,values)
    if(is.null(object$categorical)){
      output[[1L]][rows] <- as.numeric(value)
    }else if(is.matrix(output[[position]])){
      output[[position]][] <- value
    }else output[[position]] <- as.numeric(value)
  }
  if(object$singleCategorical || is.null(object$categorical)) output[[1L]] else output
}

.nonlinear_simulate <- function(object,nsim,seed,re.form,maxMemoryMB){
  nsim <- .mme_positive(nsim,"nsim",integer=TRUE)
  if(is.null(maxMemoryMB)) maxMemoryMB <- object$control$execution$denseMemoryMB
  maxMemoryMB <- .mme_positive(maxMemoryMB,"maxMemoryMB")
  if(8*as.double(object$nobs)*nsim>maxMemoryMB*1048576)
    stop("Nonlinear simulation output exceeds maxMemoryMB.",call.=FALSE)
  if(!is.null(seed)){
    if(!is.numeric(seed) || length(seed)!=1L || !is.finite(seed) || seed!=floor(seed))
      stop("seed must be NULL or one integer.",call.=FALSE)
    hadState <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
    if(hadState) oldState <- get(".Random.seed",envir=.GlobalEnv)
    on.exit(if(hadState) assign(".Random.seed",oldState,envir=.GlobalEnv) else
      if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
        rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
    set.seed(seed)
  }
  scope <- .mmes_blas_scope(object$control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  hasRandom <- ncol(object$randomDesign)>0L
  factors <- if(re.form=="population" && hasRandom)
    .laplace_random_effect_sampling_factors(object,maxMemoryMB,
      allowSingular=.nonlinear_has_reduced_rank(object$covStruct)) else NULL
  if(!isTRUE(object$nonlinearFamily)){
    object$linear.predictors <- object$fitted.values
    object$laplaceFamilies <- list(stats::gaussian())
    object$laplaceFamilyIndex <- rep(1L,length(object$input))
    object$familySpecifications <- list(list(family="gaussian",
      logParameter=log(object$dispersion)))
    object$priorWeights <- rep(1,length(object$input))
  }
  design <- .nonlinear_family_design(object,NULL,TRUE)
  simulations <- matrix(0,object$nobs,nsim)
  for(i in seq_len(nsim)){
    local <- object
    if(re.form=="population" && hasRandom){
      mode <- .laplace_draw_random_effects(factors)
      local$linear.predictors <- .nonlinear_family_eta(object,design,mode)
      local$fitted.values <- if(is.null(local$categorical))
        .nonlinear_family_output(local,local$linear.predictors) else local$fitted.values
    }
    simulations[,i] <- as.numeric(.laplace_simulate_response(local,1L,
      maxMemoryMB=maxMemoryMB))
  }
  result <- as.data.frame(simulations)
  if(isTRUE(object$singleCategorical)){
    specification <- object$familySpecifications[[1L]]
    result[] <- lapply(result,function(code){
      labels <- specification$levels
      if(specification$family=="ordinal") ordered(labels[code],levels=labels) else
        factor(labels[code],levels=labels)
    })
  }
  names(result) <- paste0("sim_",seq_len(nsim))
  result
}

.nonlinear_inference_settings <- function(object, control=NULL){
  if(is.null(control)) control <- .mmes_fitting_control(object)
  if(!inherits(control,"sommer_mmes_control"))
    stop("control must be an mmes_control() object.",call.=FALSE)
  list(numerical=if(!is.null(control$inference)) control$inference else
      .nonlinear_controls(control$task$nonlinear)$inference,
    memory=.mme_positive(control$execution$denseMemoryMB,"execution$denseMemoryMB"))
}

.nonlinear_inference_mode_settings <- function(mode){
  # Nested differences amplify mode error at perturbed outer coordinates.
  mode$scoreTol <- min(mode$scoreTol,1e-12)
  mode$stepTol <- min(mode$stepTol,1e-12)
  mode
}

.nonlinear_workspace <- function(count, memory, label){
  if(!is.finite(count) || 8*count>memory*1048576)
    stop(label," exceeds execution$denseMemoryMB.",call.=FALSE)
}

.nonlinear_block_widths <- function(object){
  if(length(object$covStruct)!=length(object$randomRelationships))
    stop("Nonlinear covariance and relationship metadata are inconsistent.",call.=FALSE)
  widths <- vapply(seq_along(object$covStruct),function(term)
    as.integer(object$covStruct[[term]]$dim*nrow(object$randomRelationships[[term]])),
    integer(1))
  if(anyNA(widths) || sum(widths)!=length(object$randomMode))
    stop("Nonlinear random-effect block metadata are inconsistent.",call.=FALSE)
  widths
}

.nonlinear_jacobian <- function(evaluate, parameters, settings, label){
  value <- as.numeric(evaluate(parameters))
  if(any(!is.finite(value))) stop(label," is non-finite.",call.=FALSE)
  calculate <- function(multiplier){
    result <- matrix(0,length(value),length(parameters))
    for(index in seq_along(parameters)){
      step <- multiplier*settings$step*max(1,abs(parameters[index]))
      upper <- lower <- parameters
      upper[index] <- upper[index]+step
      lower[index] <- lower[index]-step
      result[,index] <- (as.numeric(evaluate(upper))-as.numeric(evaluate(lower)))/(2*step)
    }
    if(any(!is.finite(result))) stop(label," is non-finite.",call.=FALSE)
    result
  }
  coarse <- calculate(1)
  fine <- calculate(0.5)
  if(any(abs(fine-coarse)>settings$derivativeTol*(1+abs(fine))))
    stop(label," did not agree at successive step sizes; revise task$nonlinear$inference$step or improve the fit.",
      call.=FALSE)
  fine
}

.nonlinear_positive_inverse <- function(information, settings, label){
  information <- as.matrix(information)
  if(any(!is.finite(information)) ||
     max(abs(information-t(information)))>
       settings$derivativeTol*(1+max(abs(information))))
    stop(label," is non-finite or asymmetric.",call.=FALSE)
  information <- (information+t(information))/2
  diagonal <- diag(information)
  if(any(!is.finite(diagonal)) || any(diagonal<=0))
    stop(label," is not positive definite; uncertainty is unavailable.",call.=FALSE)
  scale <- sqrt(diagonal)
  standardized <- information/outer(scale,scale)
  values <- eigen(standardized,symmetric=TRUE,only.values=TRUE)$values
  condition <- max(values)/min(values)
  if(min(values)<=settings$curvatureTol*max(values) ||
     !is.finite(condition) || condition>settings$maxCondition)
    stop(label," is singular, indefinite or weakly identified; uncertainty is unavailable.",
      call.=FALSE)
  result <- chol2inv(chol(standardized))/outer(scale,scale)
  if(any(!is.finite(result)))
    stop(label," inverse is outside the finite numerical range.",call.=FALSE)
  attr(result,"condition") <- condition
  result
}

.nonlinear_postfit_context <- function(object, memory){
  layout <- object$nonlinearParameterLayout$covariance
  model <- object$nonlinear
  random <- object$nonlinearRandom
  if(isTRUE(object$nonlinearGraph)){
    mode <- .nonlinear_inference_mode_settings(
      .mmes_fitting_control(object)$task$nonlinear$mode)
    if(isTRUE(object$nonlinearFamily)){
      coordinateMap <- if(is.null(random))
        matrix(0,ncol(object$randomDesign),length(object$nonlinearCoordinates)) else
          random$coordinateMap*1.0
      context <- sommer_nonlinear_graph_family_create(object$fixedDesign,
        object$randomDesign,object$nativeResponse,object$offset,object$graphData,
        object$graphNative$nodes,object$graphNative$output,
        object$graphNative$transformed,coordinateMap,object$covStruct,
        object$randomRelationships,object$familySpecifications,
        object$laplaceFamilyIndex-1L,object$priorWeights,memory,
        mode$maxIters,mode$scoreTol,mode$stepTol)
      evaluate <- function(parameters,gradient=TRUE){
        value <- .laplace_constrained_evaluate(parameters,layout,function(expanded)
          sommer_nonlinear_graph_family_eval(context,expanded,gradient))
        if(!isTRUE(value$valid))
          stop("Nonlinear graph family evaluation failed: ",value$reason,".",
            call.=FALSE)
        value
      }
      return(list(evaluate=evaluate,modeSettings=mode,
        conditional=function(parameters,latent){
        value <- sommer_nonlinear_graph_family_conditional(context,
          .laplace_expand_parameters(parameters,layout),latent)
        if(!isTRUE(value$valid))
          stop("Nonlinear graph family conditional curvature failed: ",
            value$reason,".",call.=FALSE)
        value
      }))
    }
    if(!is.null(random)){
      context <- sommer_nonlinear_graph_laplace_create(object$fixedDesign,
        object$randomDesign,object$y,object$offset,object$graphData,
        object$graphNative$nodes,object$graphNative$output,
        object$graphNative$transformed,random$coordinateMap*1.0,
        object$covStruct,object$randomRelationships,memory,
        mode$maxIters,mode$scoreTol,mode$stepTol)
      evaluate <- function(parameters,gradient=TRUE){
        value <- .laplace_constrained_evaluate(parameters,layout,function(expanded)
          sommer_nonlinear_graph_laplace_eval(context,expanded,gradient))
        if(!isTRUE(value$valid))
          stop("Nonlinear graph Laplace evaluation failed: ",value$reason,".",
            call.=FALSE)
        value
      }
      return(list(evaluate=evaluate,modeSettings=mode,
        conditional=function(parameters,latent){
        value <- sommer_nonlinear_graph_laplace_conditional(context,
          .laplace_expand_parameters(parameters,layout),latent)
        if(!isTRUE(value$valid))
          stop("Nonlinear graph conditional curvature failed: ",value$reason,".",
            call.=FALSE)
        value
      }))
    }
    context <- sommer_nonlinear_graph_gaussian_create(object$fixedDesign,
      object$randomDesign,object$y,object$offset,object$graphData,
      object$graphNative$nodes,object$graphNative$output,
      object$graphNative$transformed,object$covStruct,
      object$randomRelationships,memory)
    evaluate <- function(parameters,gradient=TRUE){
      value <- .laplace_constrained_evaluate(parameters,layout,function(expanded)
        sommer_nonlinear_gaussian_eval(context,expanded,gradient))
      if(!isTRUE(value$valid))
        stop("Nonlinear graph likelihood evaluation failed: ",value$reason,".",
          call.=FALSE)
      value
    }
    gaussian <- .laplace_family_specs(stats::gaussian(),list())$specifications
    gaussian[[1L]]$free <- FALSE
    gaussian[[1L]]$logParameter <- log(object$dispersion)
    graphMean <- .nonlinear_graph_first(model,object$trainingData,
      object$nonlinearCoordinates,object$control)$output
    mode <- .mmes_fitting_control(object)$task$nonlinear$mode
    # A zero-valued fixed curve lets the family wrapper supply Gaussian random-mode curvature.
    conditionalContext <- sommer_nonlinear_family_create(object$fixedDesign,
      object$randomDesign,object$y,object$offset+graphMean,
      rep(0,length(object$y)),"saturation",0L,ncol(object$randomDesign),
      object$covStruct,object$randomRelationships,gaussian,rep(0L,length(object$y)),
      rep(1,length(object$y)),memory,mode$maxIters,mode$scoreTol,mode$stepTol)
    p <- ncol(object$fixedDesign)
    k <- length(object$nonlinearCoordinates)
    covarianceCount <- sum(vapply(object$covStruct,function(structure)
      sum(structure$free),integer(1)))
    conditional <- function(parameters,latent){
      expanded <- .laplace_expand_parameters(parameters,layout)
      local <- c(expanded[seq_len(p)],c(0,0),
        expanded[p+k+seq_len(covarianceCount)])
      value <- sommer_nonlinear_family_conditional(conditionalContext,local,latent)
      if(!isTRUE(value$valid))
        stop("Nonlinear graph conditional curvature failed: ",value$reason,".",
          call.=FALSE)
      value
    }
    return(list(evaluate=evaluate,conditional=conditional))
  }
  mode <- .mmes_fitting_control(object)$task$nonlinear$mode
  specifications <- if(isTRUE(object$nonlinearFamily)) object$familySpecifications else
    .laplace_family_specs(stats::gaussian(),list())$specifications
  familyContext <- sommer_nonlinear_family_create(object$fixedDesign,object$randomDesign,
    if(isTRUE(object$nonlinearFamily)) object$nativeResponse else object$y,
    object$offset,object$input,model$curve,
    if(is.null(random)) 0L else random$coordinate[1L],
    if(is.null(random)) ncol(object$randomDesign) else random$start,
    object$covStruct,object$randomRelationships,specifications,
    if(isTRUE(object$nonlinearFamily)) object$laplaceFamilyIndex-1L else
      rep(0L,length(object$y)),
    if(isTRUE(object$nonlinearFamily)) object$priorWeights else rep(1,length(object$y)),
    memory,mode$maxIters,mode$scoreTol,mode$stepTol,
    coordinateMap=if(identical(random$version,2L)) random$coordinateMap*1.0 else NULL)
  context <- familyContext
  native <- sommer_nonlinear_family_eval
  if(!isTRUE(object$nonlinearFamily)){
    if(is.null(random)){
      context <- sommer_nonlinear_gaussian_create(object$fixedDesign,object$randomDesign,
        object$y,object$offset,object$input,model$curve,object$covStruct,
        object$randomRelationships,memory)
      native <- sommer_nonlinear_gaussian_eval
    }else if(!identical(random$version,2L)){
      context <- sommer_nonlinear_laplace_create(object$fixedDesign,object$randomDesign,
        object$y,object$offset,object$input,model$curve,random$coordinate,random$start,
        object$covStruct,object$randomRelationships,memory,
        mode$maxIters,mode$scoreTol,mode$stepTol)
      native <- sommer_nonlinear_laplace_eval
    }
  }
  evaluate <- function(parameters,gradient=TRUE){
    value <- .laplace_constrained_evaluate(parameters,layout,function(expanded)
      native(context,expanded,gradient))
    if(!isTRUE(value$valid))
      stop("Nonlinear post-fit likelihood evaluation failed: ",value$reason,".",call.=FALSE)
    value
  }
  list(evaluate=evaluate,conditional=function(parameters,latent){
    value <- sommer_nonlinear_family_conditional(familyContext,
      .laplace_expand_parameters(parameters,layout),latent)
    if(!isTRUE(value$valid))
      stop("Nonlinear post-fit conditional curvature failed: ",value$reason,".",call.=FALSE)
    value
  })
}

.nonlinear_at <- function(object, parameters, mode=NULL){
  p <- length(object$b)
  k <- length(object$nonlinearCoordinates)
  layout <- object$nonlinearParameterLayout$covariance
  object$b[] <- parameters[seq_len(p)]
  object$nonlinearCoordinates[] <- parameters[p+seq_len(k)]
  object$nonlinearParameters <- if(isTRUE(object$nonlinearGraph)){
    natural <- object$nonlinearCoordinates
    transformed <- which(vapply(object$nonlinear$parameterMap,function(map)
      identical(map$transform,"log"),logical(1)))
    if(length(transformed)) natural[transformed] <- exp(natural[transformed])
    natural
  }else .nonlinear_natural(object$nonlinear,object$nonlinearCoordinates)
  object$covStruct <- .laplace_proposal_structures(object$covStruct,parameters,layout)
  object$covPar <- .nonlinear_covariance_parameters(object$covStruct)
  if(isTRUE(object$nonlinearFamily)){
    object$familySpecifications <- .laplace_proposal_specs(object$familySpecifications,
      parameters,layout)
    object$laplaceFamilies <- .laplace_update_families(object$laplaceFamilies,
      object$familySpecifications,object$laplaceFamilyIndex)
  }else object$dispersion <- exp(tail(parameters,1L))
  if(!is.null(mode)){
    object$randomMode <- as.numeric(mode)
    random <- object$nonlinearRandom
    if(!is.null(random)){
      object$curveRandomMode <- mode[random$start+seq_len(random$count)]
      object$additiveRandomMode <- mode[seq_len(random$start)]
    }
  }
  object
}

.nonlinear_parameter_values <- function(object, parameters){
  local <- .nonlinear_at(object,parameters)
  covariance <- unlist(lapply(seq_along(local$covStruct),function(term){
    descriptor <- local$covStruct[[term]]
    stats::setNames(.nonlinear_covariance_parameters(list(descriptor))[[1L]],
      paste0("covariance:",local$randomMetadata[term],":",descriptor$par_names))
  }),use.names=TRUE)
  family <- c(`residual:variance`=unname(local$dispersion))
  if(isTRUE(local$nonlinearFamily)){
    family <- unlist(lapply(seq_along(local$familySpecifications),function(index){
      specification <- local$familySpecifications[[index]]
      value <- if(specification$family=="ordinal") specification$thresholds else
        if(!is.null(specification$zeroDesign)) specification$zeroCoefficients else
          if(!is.null(specification$parameter)) .laplace_family_parameter_value(specification)
          else numeric()
      if(!length(value)) return(numeric())
      stats::setNames(as.numeric(value),paste0("family:",
        names(local$familySpecifications)[index],":",
        if(specification$family=="ordinal") paste0("cutpoint",seq_along(value)) else
          if(!is.null(specification$zeroDesign)) colnames(specification$zeroDesign) else
            specification$parameter))
    }),use.names=TRUE)
  }
  c(stats::coef(local),covariance,family)
}

.nonlinear_inference_ready <- function(object){
  if(!inherits(object,"mmes.nonlinear"))
    stop("post_nonlinear() requires a nonlinear mmes fit.",call.=FALSE)
  if(!isTRUE(object$convergence) || !isTRUE(object$nonlinearDiagnostics$identified))
    stop("Nonlinear post-fit inference requires a converged, locally identifiable fit.",call.=FALSE)
  boundary <- object$nonlinearDiagnostics$boundaryParameters
  for(term in seq_along(object$covStruct)){
    if(.nonlinear_has_reduced_rank(object$covStruct[term]))
      boundary <- setdiff(boundary,paste0(term,":near-singular covariance"))
  }
  if(length(boundary))
    stop("Nonlinear uncertainty is unavailable at flagged covariance boundaries: ",
      paste(boundary,collapse=", "),". Boundary tests require calibrated bootstrap inference.",
      call.=FALSE)
}

post_nonlinear <- function(object, mode=1L, control=NULL){
  if(inherits(object,"mmes.specification"))
    return(.spec_post_nonlinear(object,mode,control))
  .nonlinear_inference_ready(object)
  if(length(mode)!=1L || is.na(mode) || !mode %in% 0:2)
    stop("mode must be 0, 1 or 2.",call.=FALSE)
  settings <- .nonlinear_inference_settings(object,control)
  numerical <- settings$numerical
  parameters <- object$optimizer$par
  d <- length(parameters); q <- length(object$randomMode)
  latent <- length(object$latentRandomMode)
  .nonlinear_workspace(20*as.double(d)^2+
    if(mode) 12*as.double(q+latent+d)^2 else 0,settings$memory,
    "Nonlinear post-fit uncertainty workspace")
  scope <- .mmes_blas_scope(object$control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  context <- .nonlinear_postfit_context(object,settings$memory)
  value <- context$evaluate(parameters)
  if(abs(value$logLik-object$llik)>1e-7*(1+abs(object$llik)))
    stop("Reconstructed nonlinear likelihood differs from the saved fit.",call.=FALSE)
  information <- -.nonlinear_jacobian(function(at) context$evaluate(at)$gradient,
    parameters,numerical,"Observed likelihood curvature")
  covariance <- .nonlinear_positive_inverse(information,numerical,"Observed likelihood curvature")
  if(max(abs(value$gradient)*sqrt(diag(covariance)))>0.1)
    stop("The nonlinear likelihood score is too large for reliable local inference; improve convergence.",
      call.=FALSE)
  dimnames(covariance) <- list(names(parameters),names(parameters))
  estimates <- .nonlinear_parameter_values(object,parameters)
  .nonlinear_workspace(20*as.double(d+length(estimates))^2+
    if(mode) 12*as.double(q+latent+d)^2 else 0,settings$memory,
    "Nonlinear natural-scale uncertainty workspace")
  transform <- .nonlinear_jacobian(function(at) .nonlinear_parameter_values(object,at),
    parameters,numerical,"Natural-scale parameter Jacobian")
  natural <- transform%*%covariance%*%t(transform)
  if(any(!is.finite(natural)) || any(diag(natural)< -1e-10))
    stop("Natural-scale parameter covariance is non-finite or has negative variance.",call.=FALSE)
  dimnames(natural) <- list(names(estimates),names(estimates))
  uncertainty <- list(version=1L,parameters=estimates,working=parameters,
    workingCovariance=covariance,parameterCovariance=natural,
    information=(information+t(information))/2,
    method=paste("Observed marginal-likelihood curvature;",
      if(isTRUE(object$nonlinearFamily) || !is.null(object$nonlinearRandom))
        "Laplace and first-order empirical Bayes delta approximations" else
          "exact Gaussian likelihood and first-order empirical Bayes delta approximation"),
    settings=settings,modeSettings=context$modeSettings,mode=as.integer(mode))
  uncertainty$outerScoreMethod <- if(!is.null(value$work$outerScore))
    as.character(value$work$outerScore) else
      if(!is.null(object$nonlinearDiagnostics$outerScore))
        as.character(object$nonlinearDiagnostics$outerScore) else "analytic"
  if(mode){
    conditional <- context$conditional(parameters,as.numeric(value$latentMode))
    inverse <- .nonlinear_positive_inverse(conditional$hessian,numerical,
      "Conditional latent-mode curvature")
    map <- sommer_nonlinear_coordinate_map(object$randomDesign,object$covStruct,
      object$randomRelationships,settings$memory)
    conditionalCovariance <- as.matrix(map%*%inverse%*%Matrix::t(map))
    motion <- .nonlinear_jacobian(function(at) context$evaluate(at)$mode,
      parameters,numerical,"Random-mode sensitivity")
    cross <- covariance%*%t(motion)
    pev <- conditionalCovariance+motion%*%cross
    if(any(!is.finite(pev)) || any(diag(pev)< -1e-10))
      stop("Nonlinear random-effect approximate PEV is non-finite or has negative variance.",
        call.=FALSE)
    uncertainty$modeSensitivity <- motion
    uncertainty$conditionalModeCovariance <- conditionalCovariance
    uncertainty$randomPEV <- if(mode==2L) pev else diag(pev)
    uncertainty$parameterModeCovariance <- if(mode==2L) cross else NULL
    uncertainty$randomMethod <- paste("Conditional covariance is inverse exact conditional curvature;",
      "PEV adds first-order estimated-parameter mode sensitivity,",
      "an empirical Bayes approximation, not an exact frequentist PEV.")
    widths <- .nonlinear_block_widths(object)
    object$randomBlockWidths <- widths
    starts <- cumsum(c(0L,head(widths,-1L)))
    object$uPevList <- lapply(seq_along(widths),function(term){
      rows <- starts[term]+seq_len(widths[term])
      levels <- nrow(object$randomRelationships[[term]])
      matrix(diag(pev)[rows],nrow=levels,
        dimnames=list(rownames(object$randomRelationships[[term]]),
          paste0(object$randomMetadata[term],".",seq_len(widths[term]/levels))))
    })
    names(object$uPevList) <- object$randomMetadata
    random <- object$nonlinearRandom
    if(!is.null(random)){
      group <- .nonlinear_jacobian(function(at){
        local <- .nonlinear_at(object,at,context$evaluate(at)$mode)
        as.numeric(.nonlinear_random_parameters(local$nonlinear,local$nonlinearCoordinates,
          random,local$randomMode))
      },parameters,numerical,"Group-parameter Jacobian")
      count <- nrow(object$groupParameters)
      groupRandom <- matrix(0,count*length(object$nonlinearCoordinates),q)
      groupMap <- if(identical(random$version,2L)) random$groupMap else{
        map <- matrix(0L,count,length(object$nonlinearCoordinates))
        map[,random$coordinate] <- random$start+seq_len(random$count)
        map
      }
      positive <- if(isTRUE(object$nonlinearGraph))
        vapply(object$nonlinear$parameterMap,function(map)
          identical(map$transform,"log"),logical(1)) else object$nonlinear$positive
      for(j in seq_len(ncol(groupMap))){
        selected <- which(groupMap[,j]>0L)
        rows <- (j-1L)*count+selected
        groupRandom[cbind(rows,groupMap[selected,j])] <-
          if(positive[j]) object$groupParameters[selected,j] else 1
      }
      groupVariance <- rowSums((group%*%covariance)*group)+
        rowSums((groupRandom%*%conditionalCovariance)*groupRandom)
      if(any(!is.finite(groupVariance)) || any(groupVariance< -1e-10))
        stop("Group-parameter delta variance is non-finite or negative.",call.=FALSE)
      uncertainty$groupParameterSE <- matrix(sqrt(pmax(0,groupVariance)),
        nrow=nrow(object$groupParameters),dimnames=dimnames(object$groupParameters))
      if(mode==2L){
        .nonlinear_workspace(12*as.double(length(object$groupParameters))^2+
          12*as.double(q+latent+d)^2,settings$memory,
          "Nonlinear group-parameter covariance workspace")
        groupCovariance <- group%*%covariance%*%t(group)+
          groupRandom%*%conditionalCovariance%*%t(groupRandom)
        labels <- as.vector(outer(rownames(object$groupParameters),
          colnames(object$groupParameters),function(group,parameter)
            paste0(parameter,"[",group,"]")))
        dimnames(groupCovariance) <- list(labels,labels)
        uncertainty$groupParameterCovariance <- groupCovariance
      }
    }
  }else object$uPevList <- NULL
  object$nonlinearInference <- uncertainty
  object$uPevMode <- as.integer(mode)
  object$nonlinearDiagnostics$uncertaintyComputed <- TRUE
  object
}

.nonlinear_uncertainty <- function(object, mode=0L){
  saved <- object$nonlinearInference
  if(is.null(saved) || saved$mode<mode)
    saved <- post_nonlinear(object,mode=mode)$nonlinearInference
  saved
}

.nonlinear_covparams_se <- function(object,term,rel_step){
  if(!identical(as.numeric(rel_step),1e-6))
    stop("Use task$nonlinear$inference$step for nonlinear reporting derivatives.",call.=FALSE)
  uncertainty <- .nonlinear_uncertainty(object)
  output <- .covparams_mmes(object,term)
  derivative <- .nonlinear_jacobian(function(at){
    reported <- .covparams_mmes(.nonlinear_at(object,at),term)
    if(!identical(reported[c("term","factor","section","parameter")],
      output[c("term","factor","section","parameter")]))
      stop("Covariance reporting changed shape under perturbation.",call.=FALSE)
    reported$estimate
  },object$optimizer$par,uncertainty$settings$numerical,"Covariance-report Jacobian")
  covariance <- derivative%*%uncertainty$workingCovariance%*%t(derivative)
  output$StdError <- sqrt(pmax(0,diag(covariance)))
  attr(output,"vcov") <- covariance
  attr(output,"method") <- uncertainty$method
  output
}

.nonlinear_covmatrix_se <- function(object,term,max.dim,rel_step){
  if(!identical(as.numeric(rel_step),1e-5))
    stop("Use task$nonlinear$inference$step for nonlinear matrix derivatives.",call.=FALSE)
  uncertainty <- .nonlinear_uncertainty(object)
  output <- covmatrix_mmes(object,term,se=FALSE,max.dim=max.dim)
  count <- length(output$covariance)
  .nonlinear_workspace(12*as.double(count)*length(object$optimizer$par),
    uncertainty$settings$memory,"Nonlinear covariance-matrix uncertainty workspace")
  derivative <- .nonlinear_jacobian(function(at){
    value <- covmatrix_mmes(.nonlinear_at(object,at),term,se=FALSE,max.dim=max.dim)
    c(as.numeric(value$covariance),as.numeric(value$correlation))
  },object$optimizer$par,uncertainty$settings$numerical,"Covariance-matrix Jacobian")
  variance <- rowSums((derivative%*%uncertainty$workingCovariance)*derivative)
  if(any(!is.finite(variance)) || any(variance< -1e-10))
    stop("Covariance-matrix delta variance is non-finite or negative.",call.=FALSE)
  output$covariance.se <- output$covariance
  output$covariance.se[] <- sqrt(pmax(0,variance[seq_len(count)]))
  output$correlation.se <- output$correlation
  output$correlation.se[] <- sqrt(pmax(0,variance[count+seq_len(count)]))
  diag(output$correlation.se) <- 0
  attr(output,"method") <- uncertainty$method
  output
}

vcDelta.mmes.nonlinear <- function(object,transform){
  if(!inherits(transform,"formula") || !length(transform) %in% c(2L,3L))
    stop("transform must be a one- or two-sided formula using V1, V2, ... .",call.=FALSE)
  uncertainty <- .nonlinear_uncertainty(object)
  selected <- which(!grepl("^(fixed|curve):",names(uncertainty$parameters)))
  values <- as.list(unname(uncertainty$parameters[selected]))
  names(values) <- paste0("V",seq_along(values))
  if(!length(values)) stop("No variance or family parameters are available.",call.=FALSE)
  value <- eval(stats::deriv(transform[[length(transform)]],names(values)),values)
  derivative <- as.numeric(attr(value,"gradient"))
  if(length(value)!=1L || !is.finite(value) ||
     length(derivative)!=length(values) || any(!is.finite(derivative)))
    stop("The variance-parameter transform or its gradient is non-finite or nonscalar.",call.=FALSE)
  covariance <- uncertainty$parameterCovariance[selected,selected,drop=FALSE]
  variance <- as.numeric(crossprod(derivative,covariance%*%derivative))
  if(!is.finite(variance) || variance< -1e-10)
    stop("The transformed variance is non-finite or negative.",call.=FALSE)
  result <- data.frame(Estimate=as.numeric(value),SE=sqrt(max(0,variance)))
  if(length(transform)==3L) rownames(result) <- paste(deparse(transform[[2L]]),collapse="")
  class(result) <- c("vcDelta.mmes","data.frame")
  attr(result,"parameters") <- stats::setNames(names(uncertainty$parameters)[selected],names(values))
  attr(result,"method") <- uncertainty$method
  result
}

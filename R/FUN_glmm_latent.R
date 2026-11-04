.laplace_latent_covariance <- function(design, structures, relationships,
                                        blockWidths, groups=NULL, memoryMB=1000){
  if(length(structures)!=length(relationships) ||
     length(structures)!=length(blockWidths) || sum(blockWidths)!=ncol(design))
    stop("Population prediction could not align latent covariance blocks.",call.=FALSE)
  if(is.null(groups)) groups <- matrix(seq_len(nrow(design)),ncol=1L)
  if(!is.matrix(groups) || !nrow(groups) || !ncol(groups) || anyNA(groups) ||
     any(groups<1 | groups>nrow(design)) || any(groups!=floor(groups)))
    stop("Population prediction latent row groups are invalid.",call.=FALSE)
  if(length(memoryMB)!=1L || !is.finite(memoryMB) || memoryMB<=0)
    stop("Latent covariance memory budget must be positive and finite.",call.=FALSE)
  d <- ncol(groups)
  bytes <- 8*as.double(nrow(groups))*d*d
  budget <- memoryMB*1048576
  if(bytes>budget)
    stop("Latent predictor covariance exceeds execution$denseMemoryMB.",call.=FALSE)
  result <- array(0,c(d,d,nrow(groups)))
  start <- 0L
  for(term in seq_along(structures)){
    descriptor <- structures[[term]]
    precision <- Matrix::Matrix(relationships[[term]],sparse=TRUE)
    coordinates <- descriptor$dim
    levels <- nrow(precision)
    width <- blockWidths[term]
    workspace <- bytes+32*as.double(coordinates)^2+24*Matrix::nnzero(precision)
    panel <- floor((budget-workspace)/(40*as.double(width)*d))
    if(!is.finite(panel) || panel<1)
      stop("Latent covariance solve workspace exceeds execution$denseMemoryMB.",call.=FALSE)
    if(width!=coordinates*levels || ncol(precision)!=levels ||
       any(!is.finite(precision)) || !isSymmetric(precision))
      stop("Population prediction latent covariance dimensions or relationship precision are invalid.",call.=FALSE)
    sigma <- evaluate_covstruct_cpp(descriptor,descriptor$par)
    if(any(!is.finite(sigma)) || !isSymmetric(sigma))
      stop("Population prediction latent coordinate covariance is invalid.",call.=FALSE)
    chol(sigma)
    precision@factors <- list()
    factor <- Matrix::Cholesky(precision,LDL=FALSE,perm=TRUE)
    block <- design[,start+seq_len(width),drop=FALSE]
    for(first in seq.int(1L,nrow(groups),by=min(panel,nrow(groups)))){
      selected <- first:min(nrow(groups),first+panel-1L)
      rows <- as.vector(t(groups[selected,,drop=FALSE]))
      rhs <- as.matrix(Matrix::t(block[rows,,drop=FALSE]))
      solved <- rhs
      for(coordinate in seq_len(coordinates)){
        positions <- (coordinate-1L)*levels+seq_len(levels)
        solved[positions,] <- as.matrix(Matrix::solve(factor,
          rhs[positions,,drop=FALSE],system="A"))
      }
      applied <- matrix(0,nrow(rhs),ncol(rhs))
      for(left in seq_len(coordinates)) for(right in seq_len(coordinates)){
        positions <- (left-1L)*levels+seq_len(levels)
        other <- (right-1L)*levels+seq_len(levels)
        applied[positions,] <- applied[positions,,drop=FALSE]+
          sigma[left,right]*solved[other,,drop=FALSE]
      }
      for(index in seq_along(selected)){
        columns <- (index-1L)*d+seq_len(d)
        value <- as.matrix(block[rows[columns],,drop=FALSE] %*%
          applied[,columns,drop=FALSE])
        result[,,selected[index]] <- result[,,selected[index]]+value
      }
    }
    start <- start+width
  }
  if(any(!is.finite(result)))
    stop("Population prediction latent covariance is non-finite.",call.=FALSE)
  result
}

.laplace_population_new_levels <- function(object, variable){
  coordinate <- any(vapply(object$covStruct,function(descriptor)
    any(vapply(descriptor$factors,function(factor){
      label <- factor$label
      !is.null(label) && variable %in% all.vars(str2lang(label))
    },logical(1))),logical(1)))
  terms <- vapply(object$randomMetadata,function(term)
    variable %in% all.vars(str2lang(term)),logical(1))
  independent <- all(vapply(object$randomRelationships[terms],function(precision)
    Matrix::isDiagonal(precision) && all(abs(diag(precision)-1)<1e-10),logical(1)))
  if(coordinate || !independent)
    stop("Population prediction at unseen latent coordinates or relationship levels is unsupported; declare the coordinate/relationship basis when fitting. New independent groups are supported.",call.=FALSE)
  invisible(NULL)
}

.laplace_population_delta <- function(object, fitted, evaluate, covariance=FALSE){
  scope <- .mmes_blas_scope(object$control$execution$blasThreads)
  on.exit(scope$restore(),add=TRUE)
  training <- object
  if(!is.null(object$categorical))
    training$familySpecifications <- object$categorical$components
  uncertainty <- .mmes_laplace_uncertainty(training)
  if(!uncertainty$valid)
    stop("Population prediction uncertainty is unavailable: ",uncertainty$reason,".",call.=FALSE)
  parameters <- object$optimizer$par
  count <- length(unlist(fitted))
  bytes <- 8*(3*as.double(count)*length(parameters)+
    if(covariance) 3*as.double(count)^2 else 0)
  if(bytes>object$laplacePostPEVSettings$denseMemoryMB*1048576)
    stop("Population prediction uncertainty exceeds execution$denseMemoryMB.",call.=FALSE)
  mean <- function(proposal){
    local <- object
    local$b <- proposal[seq_len(object$laplaceParameterLayout$fixedCount)]
    local$covStruct <- .laplace_proposal_structures(object$covStruct,proposal,
      object$laplaceParameterLayout)
    local$familySpecifications <- .laplace_proposal_specs(object$familySpecifications,
      proposal,object$laplaceParameterLayout)
    as.numeric(unlist(evaluate(local)))
  }
  jacobian <- matrix(0,count,length(parameters))
  for(position in seq_along(parameters)){
    upper <- lower <- parameters
    step <- 1e-5*max(1,abs(parameters[position]))
    upper[position] <- upper[position]+step
    lower[position] <- lower[position]-step
    jacobian[,position] <- (mean(upper)-mean(lower))/(2*step)
  }
  weighted <- jacobian %*% uncertainty$parameterCovariance
  variance <- rowSums(weighted*jacobian)
  if(any(!is.finite(variance)) || any(variance < -1e-12))
    stop("Population prediction variance is invalid; no standard error was returned.",call.=FALSE)
  standardErrors <- sqrt(pmax(variance,0))
  reshape <- function(template,values){
    if(is.matrix(template)) matrix(values,nrow(template),ncol(template),
      dimnames=dimnames(template)) else stats::setNames(values,names(template))
  }
  if(is.list(fitted)){
    starts <- cumsum(c(0L,head(lengths(fitted),-1L)))
    standardErrors <- lapply(seq_along(fitted),function(index)
      reshape(fitted[[index]],standardErrors[starts[index]+seq_along(fitted[[index]])]))
    names(standardErrors) <- names(fitted)
  }else standardErrors <- reshape(fitted,standardErrors)
  result <- list(fit=fitted,se.fit=standardErrors,uncertainty="full")
  if(covariance){
    value <- weighted %*% t(jacobian)
    result$cov.fit <- (value+t(value))/2
    labels <- .laplace_prediction_labels(fitted)
    dimnames(result$cov.fit) <- list(labels,labels)
  }
  result
}

.laplace_prediction_labels <- function(fitted){
  if(is.list(fitted)){
    return(unlist(lapply(seq_along(fitted),function(position)
      paste(names(fitted)[position],.laplace_prediction_labels(fitted[[position]]),
        sep=":")),use.names=FALSE))
  }
  if(is.matrix(fitted)) paste(rep(colnames(fitted),each=nrow(fitted)),
    rep(seq_len(nrow(fitted)),ncol(fitted)),sep=":") else as.character(seq_along(fitted))
}

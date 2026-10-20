postPEV <- function(object, mode = 1L){
  
  if(!inherits(object, "mmes")){
    stop("'object' must inherit from class 'mmes'.", call. = FALSE)
  }
  .mmes_require_coefficient_matrix(object)
  
  mode <- as.integer(mode)
  
  if(length(mode) != 1L || is.na(mode) || !mode %in% 0:2){
    stop("'mode' must be one of 0, 1, or 2.", call. = FALSE)
  }

  if(!is.null(object$rotation) && mode == 1L){
    stop(
      paste0(
        "Rotated effects require mode=2 so their full eigen-coordinate ",
        "covariance can be transformed back exactly."
      ),
      call.=FALSE
    )
  }
  
  object <- post_mme_Cinverse_cpp(object, mode)

  if(!is.null(object$rotation) && mode == 2L){
    term <- object$rotation$term
    U <- object$rotation$vectors
    ranges <- object$partitions[[term]]
    for(j in seq_len(nrow(ranges))){
      rr <- ranges[j,1]:ranges[j,2]
      block <- U %*% as.matrix(object$Ci[rr,rr,drop=FALSE]) %*% t(U)
      object$uPevList[[term]][,j] <- diag(block)
    }
  }
  
  if(length(object$uPevList) && length(object$uList)){
    names(object$uPevList) <- names(object$uList)
    
    for(i in seq_along(object$uPevList)){
      if(length(object$uPevList[[i]])){
        dimnames(object$uPevList[[i]]) <- dimnames(object$uList[[i]])
      }
    }
  }
  
  object
}


.mmes_varu_check <- function(object){
  if(!inherits(object, "mmes") || inherits(object, "mmes.glmm")){
    stop("VarU requires a Gaussian mmes fit.", call.=FALSE)
  }
  .mmes_require_coefficient_matrix(object)
  if(is.null(object$C) || !nrow(object$C) || is.null(object$Cscale)){
    stop("VarU requires a Henderson fit with stored C and Cscale; direct and task$mode='solve' fits are not supported.", call.=FALSE)
  }
}

.mmes_random_precision <- function(object){
  precision <- object$randomPrecision
  if(is.null(precision)){
    inputs <- .mmes_rebuild_inputs(object)
    if(!is.null(object$rotation) || !is.null(object$factorScoreInfo)){
      stop("Refit rotated or factor-score models to retain original relationship precisions for VarU.", call.=FALSE)
    }
    precision <- inputs$Ai[seq_along(object$uList)]
    for(index in seq_along(precision)){
      if(nrow(precision[[index]]) != nrow(object$uList[[index]])){
        stop("Rebuilt relationship precision does not match fitted effects; refit the model.", call.=FALSE)
      }
    }
  }
  precision
}

.mmes_prior_vcov <- function(object, D, precision=.mmes_random_precision(object),
                             factors=NULL){
  out <- matrix(0, nrow(D), nrow(D))
  for(index in seq_along(object$uList)){
    ranges <- object$partitions[[index]]
    columns <- unlist(lapply(seq_len(nrow(ranges)), function(coordinate){
      seq.int(ranges[coordinate,1], ranges[coordinate,2])
    }), use.names=FALSE)
    if(!any(D[,columns,drop=FALSE] != 0)) next
    Sigma <- object$theta[[index]]
    dimension <- object$covStruct[[index]]$dim
    if(!identical(dim(Sigma), c(dimension, dimension))){
      Sigma <- covmatrix_mmes(object, index, se=FALSE, max.dim=dimension)$covariance
    }
    blocks <- lapply(seq_len(nrow(ranges)), function(coordinate){
      D[,seq.int(ranges[coordinate,1], ranges[coordinate,2]),drop=FALSE]
    })
    factor <- if(is.null(factors)){
      Matrix::Cholesky(Matrix::forceSymmetric(precision[[index]]), LDL=FALSE)
    }else factors[[index]]
    solved <- Matrix::solve(factor, do.call(cbind, lapply(blocks, t)), system="A")
    for(first in seq_along(blocks)){
      for(second in seq_along(blocks)){
        rhsColumns <- (second-1L)*nrow(D) + seq_len(nrow(D))
        out <- out + Sigma[first,second] * as.matrix(blocks[[first]] %*% solved[,rhsColumns,drop=FALSE])
      }
    }
  }
  (out+t(out))/2
}

.mmes_vcov_solver <- function(object){
  .mmes_require_coefficient_matrix(object)
  factor <- mmes_vcov_factor_cpp(object)
  function(D){
    mmes_vcov_solve_cpp(factor, to_sparse(.mmes_engine_contrast(object, D)))
  }
}

.mmes_henderson_vcov <- function(object, contrasts){
  solve <- .mmes_vcov_solver(object)
  lapply(contrasts, solve)
}

postVarU <- function(object, mode=1L){
  .mmes_varu_check(object)
  if(length(mode) != 1L || is.na(mode) || !mode %in% 0:2){
    stop("mode must be 0, 1, or 2.", call.=FALSE)
  }
  if(mode == 0L){
    object$uVarList <- NULL
    object$VarU <- NULL
    object$uVarMode <- 0L
    return(object)
  }
  precision <- .mmes_random_precision(object)
  priorFactors <- lapply(precision, function(x)
    Matrix::Cholesky(Matrix::forceSymmetric(x), LDL=FALSE))
  solve <- .mmes_vcov_solver(object)
  nFixed <- nrow(object$b)
  nRandom <- nrow(object$bu)-nFixed
  diagonal <- numeric(nRandom)
  if(mode == 2L){
    D <- cbind(Matrix::Matrix(0, nRandom, nFixed, sparse=TRUE), Matrix::Diagonal(nRandom))
    error <- solve(D)
    object$VarU <- .mmes_prior_vcov(object, D, precision, priorFactors)-error
    labels <- unlist(lapply(seq_along(object$uList), function(index){
      effects <- object$uList[[index]]
      paste(names(object$uList)[index],
            rep(colnames(effects), each=nrow(effects)),
            rep(rownames(effects), times=ncol(effects)), sep=":")
    }))
    dimnames(object$VarU) <- list(labels, labels)
    diagonal <- diag(object$VarU)
  }else{
    object$VarU <- NULL
    for(start in seq.int(1L, max(nRandom,1L), by=128L)){
      if(!nRandom) break
      rows <- seq.int(start, min(start+127L, nRandom))
      D <- Matrix::sparseMatrix(i=seq_along(rows), j=nFixed+rows, x=1,
                                dims=c(length(rows), nrow(object$bu)))
      error <- solve(D)
      diagonal[rows] <- diag(.mmes_prior_vcov(object, D, precision, priorFactors)-error)
    }
  }
  object$uVarList <- lapply(seq_along(object$uList), function(index){
    ranges <- object$partitions[[index]]
    columns <- unlist(lapply(seq_len(nrow(ranges)), function(coordinate){
      seq.int(ranges[coordinate,1], ranges[coordinate,2])
    }), use.names=FALSE)-nFixed
    matrix(diagonal[columns], nrow=nrow(object$uList[[index]]),
           dimnames=dimnames(object$uList[[index]]))
  })
  names(object$uVarList) <- names(object$uList)
  object$uVarMode <- as.integer(mode)
  object
}

corImputation <- function(wide, Gu=NULL, nearest=10, roundR=FALSE){
  if(is.null(rownames(wide))){stop("Rownames of the input matrix cannot be NULL. Please add them", call. = FALSE)}
  if(is.null(Gu)){
    X <- apply(wide, 2, imputev)
    Gu <- cor(t(X))
  }
  wide2 <- wide
  rowNamesWide <-  rownames(wide)
  # for each feature
  for(iEnv in 1:ncol(wide)){ # iEnv=10
    withData <- which(!is.na(wide[,iEnv]))
    withoutData <- which(is.na(wide[,iEnv]))
    toPredict <- 1:nrow(wide)
    # if(length(toPredict) > 1){
    likelihood=Gu[as.character(rowNamesWide)[toPredict],as.character(rowNamesWide)[withData], drop=FALSE]
    # }else{
    #   likelihood=matrix(Gu[as.character(rowNamesWide)[toPredict],as.character(rowNamesWide)[withData]], nrow=1)
    #   rownames(likelihood) <- as.character(rowNamesWide)[toPredict]
    #   colnames(likelihood) <- as.character(rowNamesWide)[withData]
    # }
    replacement <- numeric()
    for(iInd in 1:nrow(likelihood)){ # iInd=1
      # averaging only the positively correlated to avoid decrease   # wide[iInd, iEnv]
      indLik <- sort(abs(likelihood[iInd,]), decreasing = TRUE)
      toAverage <- indLik[1:min(c(nearest,length(withData)))]
      indLikToAverage <- likelihood[iInd,names(toAverage)]
      replacement[iInd] <- mean(wide[names(which(indLikToAverage > 0)),iEnv]) 
    }
    names(replacement) <- rownames(likelihood)
    # time to replace the missing data
    dd=data.frame(replacement=replacement, index=1:length(replacement),
                  imputed=1, id=rownames(likelihood),orVal=wide[,iEnv])
    dd$imputed[which(dd$id %in% names(withoutData))]=0
    head(dd)
    if(roundR){
      wide2[toPredict,iEnv] <- round(replacement)
    }else{
      wide2[toPredict,iEnv] <- replacement
    }
  }
  # for each individual
  for(jRow in 1:nrow(wide)){
    miss <- which(is.na(wide[jRow,]))
    if(length(miss) > 0){
      dd=data.frame(full=as.vector(unlist(wide2[jRow,])), partial=as.vector(unlist(wide[jRow,])))
      # model <- fastLm(partial~full,data=dd[which(!is.na(dd$partial)),])
      # model <- lm(partial~full,data=dd[which(!is.na(dd$partial)),])
      model <- mmes(partial~full,data=dd[which(!is.na(dd$partial)),],
                    control=mmes_control(diagnostics=list(verbose=FALSE)))
      pp=as.vector(model$b[1,1])+(dd[which(is.na(dd$partial)),"full"]*as.vector(model$b[2,1]))
      if(roundR){
        wide[jRow,miss] <- round(pp)#round(predict(model,newdata = dd[which(is.na(dd$partial)),]))
      }else{
        wide[jRow,miss] <- pp#predict(model,newdata = dd[which(is.na(dd$partial)),])
      }
    }
  }
  stillEmpty <- which(is.na(wide), arr.ind = TRUE)
  if(nrow(stillEmpty) > 0){wide[stillEmpty] <- mean(wide, na.rm=TRUE)}
  stillEmpty <- which(is.na(wide2), arr.ind = TRUE)
  if(nrow(stillEmpty) > 0){wide2[stillEmpty] <- mean(wide2, na.rm=TRUE)}
  
  return(list(imputed=wide, corImputed=wide2))
}


r2 <- function(object, object2=NULL){
  if(!inherits(object, "mmes")){
    stop("This function is only available for models fitted with the mmes() function.", call. = FALSE)
  }
  result <- list()
  for(iPart in 1:length(object$uPevList)){
    pev <- object$uPevList[[iPart]]
    variance0 <- as.numeric(diag(object$theta[[iPart]]))
    variance <- apply(data.frame(variance0),1,function(x){rep(x,nrow(pev))})
    if(!is.null(object2)){
      if(!is.null(object2$Ai)){
        relationshipDiagonal <- diag(solve(object2$Ai[[iPart]]))
        for(iVar in 1:length(variance0)){ # if user provided A matrices we subsitute with more accurate values
          variance[,iVar] <- variance0[iVar] * relationshipDiagonal
        }
      }else{
        stop("object2 needs to be a model that sets the argument 'returnParam=TRUE' so we can extract the relationship matrices. Please correct. ", call. = FALSE)
      }
      
    }
    result[[iPart]] <- (variance - pev)/variance
  }
  names(result) <- names(object$uList)
  return(result)
}




matrix.trace <- function(x){
  if (!is.square.matrix(x))
    stop("argument x is not a square matrix")
  return(sum(diag(x)))
}

is.diagonal.matrix <- function (x, tol = 1e-08){
  y <- x
  diag(y) <- rep(0, nrow(y))
  return(all(abs(y) < tol))
}

is.square.matrix <-function(x){
  return(nrow(x) == ncol(x))
}

hadamard.prod <-function (x, y){
  if (!is.numeric(x)) {
    stop("argument x is not numeric")
  }
  if (!is.numeric(y)) {
    stop("argument y is not numeric")
  }
  if (is.matrix(x)) {
    Xmat <- x
  }
  else {
    if (is.vector(x)) {
      Xmat <- matrix(x, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (is.matrix(y)) {
    Ymat <- y
  }
  else {
    if (is.vector(y)) {
      Ymat <- matrix(y, nrow = length(x), ncol = 1)
    }
    else {
      stop("argument x is neither a matrix or a vector")
    }
  }
  if (nrow(Xmat) != nrow(Ymat))
    stop("argumentx x and y do not have the same row order")
  if (ncol(Xmat) != ncol(Ymat))
    stop("arguments x and y do not have the same column order")
  return(Xmat * Ymat)
}





# Least squares for variance-component starting values only. A well-conditioned
# full-rank design goes through the normal equations (sparse crossprod plus a
# Cholesky factor); anything else uses the dense QR, whose NA coefficients flag
# aliased columns to the caller.
.mmes_start_lsq <- function(X, Y){
  p <- ncol(X)
  XtX <- as.matrix(Matrix::crossprod(X))
  R <- tryCatch(chol(XtX), error=function(e) NULL)
  if(!is.null(R) && all(is.finite(R)) &&
     min(diag(R) / sqrt(diag(XtX))) > 1e-7){
    XtY <- as.matrix(Matrix::crossprod(X, Y))
    beta <- backsolve(R, backsolve(R, XtY, transpose=TRUE))
    return(list(beta=beta, rank=p, resid=Y - as.matrix(X %*% beta)))
  }
  Xd <- as.matrix(X)
  qrX <- qr(Xd)
  beta <- if(nrow(Xd) - qrX$rank >= 1L) qr.coef(qrX, Y) else NULL
  list(beta=beta, rank=qrX$rank,
       resid=if(is.null(beta) || anyNA(beta)) NULL else Y - Xd %*% beta)
}

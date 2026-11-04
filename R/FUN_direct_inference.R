.mmes_direct_memory <- function(object, maxMemoryMB=NULL){
  if(is.null(maxMemoryMB)) maxMemoryMB <- object$control$execution$denseMemoryMB
  if(is.null(maxMemoryMB)) maxMemoryMB <- object$engineControl$denseMemoryMB
  if(is.null(maxMemoryMB)) maxMemoryMB <- 1000
  .mme_positive(maxMemoryMB, "maxMemoryMB")
}

.mmes_require_direct_inference <- function(object){
  .mmes_require_coefficient_matrix(object)
  if(is.null(object$directInference)){
    stop("Direct post-fit covariance inputs are unavailable; refit with output=list(coefficientMatrix=TRUE).",
      call.=FALSE)
  }
  invisible(NULL)
}

.mmes_direct_vcov_solver <- function(object, maxMemoryMB=NULL){
  .mmes_require_direct_inference(object)
  factor <- direct_mmes_vcov_factor_cpp(object$directInference,
    .mmes_direct_memory(object, maxMemoryMB))
  function(D, prior=FALSE){
    if(ncol(D) != nrow(object$bu)){
      stop("Prediction contrasts must have one column per fitted coefficient.", call.=FALSE)
    }
    labels <- rownames(D)
    if(!is.null(object$coefficientMap)) D <- D %*% object$coefficientMap
    result <- direct_mmes_vcov_solve_cpp(factor, to_sparse(D), prior)
    if(prior) dimnames(result) <- list(labels, labels)
    result
  }
}

.mmes_direct_postPEV <- function(object, mode){
  count <- nrow(object$bu)
  if(mode == 0L){
    object$Ci <- Matrix::sparseMatrix(i=integer(), j=integer(), dims=c(0L, 0L))
    object$uPevList <- lapply(object$uList, function(x) matrix(numeric(), 0L, 0L))
  }else{
    budget <- .mmes_direct_memory(object)
    if(mode == 2L && 8e-6 * 8 * count^2 > budget){
      stop("Full direct coefficient covariance exceeds execution$denseMemoryMB; use postPEV(mode=1) or increase the memory budget.",
        call.=FALSE)
    }
    solve <- .mmes_direct_vcov_solver(object, budget)
    if(mode == 2L){
      object$Ci <- solve(Matrix::Diagonal(count))
      diagonal <- diag(object$Ci)
    }else{
      diagonal <- numeric(count)
      panelSize <- object$control$output$panelSize
      if(is.null(panelSize)) panelSize <- 64L
      for(first in seq.int(1L, count, by=panelSize)){
        rows <- seq.int(first, min(first + panelSize - 1L, count))
        selection <- Matrix::sparseMatrix(i=seq_along(rows), j=rows, x=1,
          dims=c(length(rows), count))
        diagonal[rows] <- diag(solve(selection))
      }
      object$Ci <- Matrix::Diagonal(count, x=diagonal)
    }
    object$uPevList <- lapply(seq_along(object$uList), function(term){
      ranges <- object$partitions[[term]]
      columns <- unlist(lapply(seq_len(nrow(ranges)), function(coordinate)
        seq.int(ranges[coordinate,1L], ranges[coordinate,2L])), use.names=FALSE)
      matrix(diagonal[columns], nrow=nrow(object$uList[[term]]),
        dimnames=dimnames(object$uList[[term]]))
    })
  }
  names(object$uPevList) <- names(object$uList)
  object$CiMode <- as.integer(mode)
  object$CiComputed <- mode == 2L
  object
}

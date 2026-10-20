.relationship_impute <- function(X){
  if(!is.matrix(X) || !is.numeric(X) || nrow(X) < 2L || ncol(X) < 1L){
    stop("X must be a numeric matrix with at least two individuals and one marker.", call.=FALSE)
  }
  missingCounts <- colSums(is.na(X))
  uniqueCols <- which(missingCounts > 0L)
  if(length(uniqueCols)){
    if(any(missingCounts == nrow(X))){
      stop("Cannot mean-impute a marker with no observed genotypes.", call.=FALSE)
    }
    cat("Imputing markers with mean value\n")
    for(j in uniqueCols){
      missing <- is.na(X[,j])
      X[missing,j] <- mean(X[,j], na.rm=TRUE)
    }
  }
  X
}

A.mat <- function(X, min.MAF=0, return.imputed=FALSE){
  X <- .relationship_impute(X)
  vanraden=TRUE
  ##################
  res <- .Call("_sommer_amat",PACKAGE = "sommer",
               X, vanraden, min.MAF
  )
  colnames(res) <- rownames(res) <- rownames(X)
  if(return.imputed){
    return(list(X=X,A=res))
  }else{
    return(res)
  }
  
}

D.mat <- function(X, nishio=TRUE, min.MAF=0, return.imputed=FALSE){
  X <- .relationship_impute(X)
  ##################
  res <- .Call("_sommer_dmat",PACKAGE = "sommer", 
               X, nishio, min.MAF
  )
  colnames(res) <- rownames(res) <- rownames(X)
  if(return.imputed){
    return(list(X=X,D=res))
  }else{
    return(res)
  }
}

APY <- function(G, core, tol=1e-10, return.details=FALSE){
  G <- as.matrix(G)
  if(!is.numeric(G) || length(dim(G)) != 2L || nrow(G) != ncol(G) ||
     nrow(G) < 1L || any(!is.finite(G))){
    stop("G must be a finite, non-empty square numeric relationship matrix.", call.=FALSE)
  }
  if(length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1){
    stop("tol must be one finite value strictly between zero and one.", call.=FALSE)
  }
  scale <- max(1, max(abs(diag(G))))
  if(max(abs(G - t(G))) > tol * scale){
    stop("G must be symmetric within the requested tolerance.", call.=FALSE)
  }
  G <- 0.5 * (G + t(G))

  ids <- rownames(G)
  if(!is.null(ids) && (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))){
    stop("G row names must be non-missing, non-empty and unique.", call.=FALSE)
  }
  if(!is.null(colnames(G)) && !identical(colnames(G), ids)){
    stop("G row and column names must be identical and in the same order.", call.=FALSE)
  }
  if(is.character(core)){
    if(is.null(ids) || anyNA(core) || anyDuplicated(core) || any(!core %in% ids)){
      stop("Character core identifiers must be unique row names of G.", call.=FALSE)
    }
    coreIndex <- match(core, ids)
  }else if(is.numeric(core) || is.integer(core)){
    if(length(core) < 1L || any(!is.finite(core)) || any(core != as.integer(core)) ||
       any(core < 1L | core > nrow(G)) || anyDuplicated(core)){
      stop("core must contain unique valid one-based row indices or row names.", call.=FALSE)
    }
    coreIndex <- as.integer(core)
  }else{
    stop("core must contain row names or one-based row indices.", call.=FALSE)
  }

  nonCoreIndex <- setdiff(seq_len(nrow(G)), coreIndex)
  Gcc <- G[coreIndex, coreIndex, drop=FALSE]
  coreChol <- tryCatch(chol(Gcc), error=function(e) NULL)
  if(is.null(coreChol)){
    stop("The core relationship submatrix is not positive definite.", call.=FALSE)
  }
  GccInverse <- chol2inv(coreChol)
  if(length(nonCoreIndex)){
    Gnc <- G[nonCoreIndex, coreIndex, drop=FALSE]
    B <- t(backsolve(coreChol, forwardsolve(t(coreChol), t(Gnc))))
    conditionalVariance <- diag(G)[nonCoreIndex] - rowSums(B * Gnc)
    if(any(!is.finite(conditionalVariance)) ||
       any(conditionalVariance <= tol * scale)){
      stop("APY conditional variances must exceed tol times the relationship scale; revise the core or explicitly blend G before APY().", call.=FALSE)
    }
    weightedB <- B / conditionalVariance
    Qcc <- GccInverse + crossprod(B, weightedB)
    Qcn <- -t(weightedB)
  }else{
    conditionalVariance <- numeric()
    B <- matrix(numeric(), 0L, length(coreIndex))
    Qcc <- GccInverse
    Qcn <- matrix(numeric(), length(coreIndex), 0L)
  }

  crossRows <- rep(coreIndex, times=length(nonCoreIndex))
  crossColumns <- rep(nonCoreIndex, each=length(coreIndex))
  Q <- Matrix::drop0(Matrix::sparseMatrix(
    i=c(rep(coreIndex, times=length(coreIndex)),
        crossRows, crossColumns, nonCoreIndex),
    j=c(rep(coreIndex, each=length(coreIndex)),
        crossColumns, crossRows, nonCoreIndex),
    x=c(as.vector(Qcc), as.vector(Qcn), as.vector(Qcn),
        1 / conditionalVariance),
    dims=dim(G), dimnames=list(ids, ids)
  ))
  attr(Q, "inverse") <- TRUE
  attr(Q, "APY") <- list(
    core=if(is.null(ids)) coreIndex else ids[coreIndex],
    coreIndex=coreIndex,
    nonCore=if(is.null(ids)) nonCoreIndex else ids[nonCoreIndex],
    conditionalVariance=conditionalVariance,
    tolerance=tol
  )
  if(return.details){
    return(list(Gu=Q, core=coreIndex, nonCore=nonCoreIndex,
                B=B, conditionalVariance=conditionalVariance))
  }
  Q
}

E.mat <- function(X,nishio=TRUE,type="A#A",min.MAF=0.02){
  if(length(type) != 1L || is.na(type) || !type %in% c("A#A", "A#D", "D#D")){
    stop("type must be one of 'A#A', 'A#D' or 'D#D'.", call.=FALSE)
  }
  X <- .relationship_impute(X)
  E <- .Call("_sommer_emat_genotypes",PACKAGE = "sommer",
             X, nishio, type, min.MAF)
  colnames(E) <- rownames(E) <- rownames(X)
  return(E)
}

H.mat <- function(A, G, tau = 1, omega = 1, tolparinv=1e-6){
  for(name in c("A", "G")){
    M <- if(name == "A") A else G
    ids <- rownames(M)
    if(!is.matrix(M) || !is.numeric(M) || nrow(M) < 1L || nrow(M) != ncol(M)){
      stop(name, " must be a non-empty square numeric matrix.", call.=FALSE)
    }
    if(is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
       !identical(colnames(M), ids)){
      stop(name, " must have identical, non-missing, non-empty and unique row and column names.", call.=FALSE)
    }
  }
  idA <- rownames(A)
  idG <- rownames(G)
  if(!all(idG %in% idA)){
    stop("All individuals in G must also occur in A.", call.=FALSE)
  }
  idH <- unique(c(idG, idA))
  idH <- rev(idH)
  A <- A[idH, idH, drop=FALSE]
  index = is.na(match(idH, idG))
  G22 <- G[idH[!index], idH[!index], drop=FALSE]
  H <- .Call("_sommer_hmat",PACKAGE = "sommer",
             A, G22, index, tolparinv, tau, omega
  )
  colnames(H) <- rownames(H) <- idH
  return(H)
  
}

dfToMatrix <- function(x, row="Row",column="Column",value="Ainverse", returnInverse=FALSE, bend=1e-6){
  
  # if(type=="asreml"){
  all <- unique(c(as.character(x[,row]),as.character(x[,column])))
  alldf <- data.frame(x=1:length(all)); rownames(alldf) <- all
  x[,row] <- as.numeric(as.character(alldf[x[,row],"x"]))
  x[,column] <- as.numeric(as.character(alldf[x[,column],"x"]))
  
  K <- matrix(NA,max(x[,row]),max(x[,column]))
  for(i in 1:nrow(x)){
    K[x[i,row],x[i,column]] <- x[i,value]
  }
  # K[1:3,1:3]
  copying <- function(m) { # copy upper triangular in lower triangular
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  copying2 <- function(m) { # copy lower triangular in upper triangular
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    m
  }
  
  K <- copying2(K)
  K[which(is.na(K), arr.ind = TRUE)] <- 0
  
  rownames(K) <- colnames(K) <- attr(x, "rowNames")
  
  Ks <- as(as(as( K,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
  if(returnInverse){
    Ksi <- solve(Ks + diag(bend, nrow(Ks)))
    rownames(Ksi) <- colnames(Ksi) <- attr(x, "rowNames")
  }else{
    Ksi <- NULL
  }
  return(list(K=Ks, Kinv=Ksi))
  # }
}

simGECorMat <- function(nEnv,nMegaEnv, mu=0.7, v=0.2, mu2=0, v2=0.3){
  ff <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  
  G = matrix(NA,nEnv,nEnv)
  (nEnv2 <- nEnv/nMegaEnv)
  G
  
  starts <- seq(1,nEnv,nEnv/nMegaEnv)
  ends <- c((starts-1)[-1],nEnv)
  
  for(i in 1:nMegaEnv){
    corsprov <- rnorm((nEnv2*(nEnv2-1))/2,mu,v)
    counter=1
    for(j in starts[i]:ends[i]){
      for(k in j:ends[i]){
        if(j == k){
          G[j,k] <- 1
        }else{
          G[j,k] <- corsprov[counter]
          counter <- counter+1
        }
      }
    }
  }
  G <- ff(G)
  tofill <- which(is.na(G),arr.ind = TRUE)
  G[tofill] <- rnorm(nrow(tofill),mu2,v2)
  G[which((G) > 1)] <- .98
  G[which((G) < -1)] <- -.98
  G <- ff(G)
  return(G)
}


AR1 = function(x,rho=0.25) {
  unx <- levels(as.factor(x))
  dims <- length(unx)
  M = diag(dims)
  M = rho^abs(row(M)-col(M))
  colnames(M) <- rownames(M) <- unx
  return(M)
}

CS = function(x, rho=0.25) {
  unx <- levels(as.factor(x))
  dims <- length(unx)
  M = matrix(rho,dims,dims)
  diag(M) <- 1
  colnames(M) <- rownames(M) <- unx
  return(M)
}

ARMA = function(x, rho=0.25, lambda=0.25) {
  ## for ar
  unx <- levels(as.factor(x))
  dimo <- length(unx)
  M = diag(dimo)
  M = abs(row(M)-col(M))
  M[lower.tri(M)] <- M[lower.tri(M)]-1
  M[upper.tri(M)] <- M[upper.tri(M)]-1
  MM <- rho^M
  ## for lam
  N <- matrix(lambda,dimo,dimo)
  diag(N) <- 1 # or 0?
  ## final
  MN <- MM*N
  colnames(MN) <- rownames(MN) <- unx
  return(MN)
}

A.mat <- function(X, endelman=TRUE, min.MAF=0, return.imputed=FALSE){
  ################
  ## impute
  missingCheck <- which(is.na(X), arr.ind = TRUE)
  if(nrow(missingCheck) > 0){
    cat("Imputing markers with mean value\n")
    uniqueCols <- unique(missingCheck[,2])
    X[,uniqueCols] <- apply(X[,uniqueCols],2,imputev)
  }
  
  ##################
  res <- .Call("_sommer_amat",PACKAGE = "sommer",
               X, endelman, min.MAF
  )
  colnames(res) <- rownames(res) <- rownames(X)
  if(return.imputed){
    return(list(X=X,A=res))
  }else{
    return(res)
  }
  
}

D.mat <- function(X, nishio=TRUE, min.MAF=0, return.imputed=FALSE){
  ################
  ## impute
  missingCheck <- which(is.na(X), arr.ind = TRUE)
  if(nrow(missingCheck) > 0){
    cat("Imputing markers with mean value\n")
    uniqueCols <- unique(missingCheck[,2])
    X[,uniqueCols] <- apply(X[,uniqueCols],2,imputev)
  }
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

E.mat <- function(X,endelman=TRUE,nishio=TRUE,type="A#A",min.MAF=0.02){
  
  if(type == "A#A"){
    A <- A.mat(X, endelman=endelman,min.MAF=min.MAF)
    E <- .Call("_sommer_emat",PACKAGE = "sommer",
               A, A
    )
  }
  if(type == "A#D"){
    A <- A.mat(X, endelman=endelman,min.MAF=min.MAF)
    D <- D.mat(X, nishio=nishio,min.MAF=min.MAF)
    E <- .Call("_sommer_emat",PACKAGE = "sommer",
               A, D
    )
  }
  if(type == "D#D"){
    D <- D.mat(X, nishio=nishio,min.MAF=min.MAF)
    E <- .Call("_sommer_emat",PACKAGE = "sommer",
               D, D
    )
  }
  colnames(E) <- rownames(E) <- rownames(X)
  return(E)
}

H.mat <- function(A, G, tau = 1, omega = 1, tolparinv=1e-6){
  
  idA <- rownames(A)
  idG <- rownames(G)
  idH <- unique(c(idG, idA))
  idH <- rev(idH)
  A <- A[idH, idH]
  index = is.na(match(idH, idG))
  G22 <- G[idH[!index], idH[!index]]
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
  
  Ks <- as(K, Class = "dgCMatrix")
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

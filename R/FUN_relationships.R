# AR1.mat = function(rho, dims) {
#   M = diag(dims)
#   M = rho^abs(row(M)-col(M))
#   return(M)
# }

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

# ARMA.mat = function(ar, lam, dimo) {
#   ## for ar
#   M = diag(dimo)
#   M = abs(row(M)-col(M))
#   M[lower.tri(M)] <- M[lower.tri(M)]-1
#   M[upper.tri(M)] <- M[upper.tri(M)]-1
#   MM <- ar^M
#   ## for lam
#   N <- matrix(lam,dimo,dimo)
#   diag(N) <- 0
#   ## final
#   MN <- MM*N
#   return(MN)
# }

A.mat <- function(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,n.core=1,shrink=FALSE,return.imputed=FALSE, ploidy=2){
  if(ploidy == 2){
    impute.EM <- function(W, cov.mat, mean.vec) {
      n <- nrow(W)
      m <- ncol(W)
      S <- matrix(0,n,n)
      for (i in 1:m) {
        Wi <- matrix(W[,i],n,1)
        missing <- which(is.na(Wi))
        if (length(missing) > 0) {
          not.NA <- setdiff(1:n,missing)
          Bt <- solve(cov.mat[not.NA,not.NA],cov.mat[not.NA,missing])
          Wi[missing] <- mean.vec[missing] + crossprod(Bt,Wi[not.NA]-mean.vec[not.NA])
          C <- cov.mat[missing,missing] - crossprod(cov.mat[not.NA,missing],Bt)
          D <- tcrossprod(Wi)
          D[missing,missing] <- D[missing,missing] + C
          W[,i] <- Wi
        } else {D <- tcrossprod(Wi)}
        S <- S + D
      }	
      return(list(S=S,W.imp=W))
    }
    
    cov.W.shrink <- function(W) {
      m <- ncol(W)
      n <- nrow(W)
      Z <- t(scale(t(W),scale=FALSE))
      Z2 <- Z^2
      S <- tcrossprod(Z)/m
      target <- mean(diag(S))*diag(n)
      var.S <- tcrossprod(Z2)/m^2-S^2/m
      b2 <- sum(var.S)
      d2 <- sum((S-target)^2)
      delta <- max(0,min(1,b2/d2))
      print(paste("Shrinkage intensity:",round(delta,2)))
      return(target*delta + (1-delta)*S)
    }
    
    X <- as.matrix(X)
    n <- nrow(X)
    frac.missing <- apply(X,2,function(x){length(which(is.na(x)))/n})
    missing <- max(frac.missing) > 0
    freq <- apply(X + 1, 2, function(x) {mean(x, na.rm = missing)})/2
    MAF <- apply(rbind(freq,1-freq),2,min)
    if (is.null(min.MAF)) {min.MAF <- 1/(2*n)}
    if (is.null(max.missing)) {max.missing <- 1 - 1/(2*n)}
    markers <- which((MAF >= min.MAF)&(frac.missing <= max.missing)) 
    m <- length(markers)
    var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
    one <- matrix(1, n, 1)
    
    mono <- which(freq*(1-freq)==0)
    X[,mono] <- 2*tcrossprod(one,matrix(freq[mono],length(mono),1))-1
    
    freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
    W <- X[, markers] + 1 - 2 *freq.mat 
    
    if (!missing) {
      if (shrink) {
        W.mean <- rowMeans(W)
        cov.W <- cov.W.shrink(W)
        A <- (cov.W+tcrossprod(W.mean))/var.A	
      } else {
        A <- tcrossprod(W)/var.A/m	
      }
      rownames(A) <- rownames(X)
      colnames(A) <- rownames(A)
      if (return.imputed) {
        return(list(A=A,imputed=X))		
      } else {
        return(A)
      }
    } else {
      #impute
      isna <- which(is.na(W)) 
      W[isna] <- 0
      
      if (toupper(impute.method)=="EM") {
        if (m < n) {
          print("Linear dependency among the lines: imputing with mean instead of EM algorithm.")
        } else {
          mean.vec.new <- matrix(rowMeans(W),n,1)
          cov.mat.new <- cov(t(W))
          if (qr(cov.mat.new)$rank < nrow(cov.mat.new)-1) {
            print("Linear dependency among the lines: imputing with mean instead of EM algorithm.")
          } else {
            
            #do EM algorithm
            W[isna] <- NA
            A.new <- (cov.mat.new + tcrossprod(mean.vec.new))/var.A
            err <- tol+1
            print("A.mat converging:")
            while (err >= tol) {
              A.old <- A.new
              cov.mat.old <- cov.mat.new
              mean.vec.old <- mean.vec.new
              if ((n.core > 1) & requireNamespace("parallel",quietly=TRUE)) {
                it <- split(1:m,factor(cut(1:m,n.core,labels=FALSE)))
                pieces <- parallel::mclapply(it,function(mark2){impute.EM(W[,mark2],cov.mat.old,mean.vec.old)},mc.cores=n.core)
              } else {
                pieces <- list()
                pieces[[1]] <- impute.EM(W,cov.mat.old,mean.vec.old)
              }
              n.pieces <- length(pieces)
              S <- matrix(0,n,n)
              W.imp <- numeric(0)
              for (i in 1:n.pieces) {
                S <- S + pieces[[i]]$S
                W.imp <- cbind(W.imp,pieces[[i]]$W.imp)
              }
              mean.vec.new <- matrix(rowMeans(W.imp),n,1)
              cov.mat.new <- (S-tcrossprod(mean.vec.new)*m)/(m-1)
              A.new <- (cov.mat.new + tcrossprod(mean.vec.new))/var.A
              err <- norm(A.old-A.new,type="F")/n
              print(err,digits=3)
            }
            rownames(A.new) <- rownames(X)
            colnames(A.new) <- rownames(A.new)
            
            if (return.imputed) {
              Ximp <- W.imp - 1 + 2*freq.mat
              colnames(Ximp) <- colnames(X)[markers]
              rownames(Ximp) <- rownames(X)
              return(list(A=A.new,imputed=Ximp))
            } else {
              return(A.new)
            }
          } #else EM 
        } #else EM
      } #else EM
      
      #imputing with mean
      if (shrink) {
        W.mean <- rowMeans(W)
        cov.W <- cov.W.shrink(W)
        A <- (cov.W+tcrossprod(W.mean))/var.A	
      } else {
        A <- tcrossprod(W)/var.A/m	
      }
      rownames(A) <- rownames(X)
      colnames(A) <- rownames(A)
      
      if (return.imputed) {
        Ximp <- W - 1 + 2*freq.mat
        colnames(Ximp) <- colnames(X)[markers]
        rownames(Ximp) <- rownames(X)
        return(list(A=A,imputed=Ximp))		
      } else {
        return(A)
      }
    } #else missing
    
  }else{
    M <- scale(X, center = TRUE, scale = FALSE)
    K <- tcrossprod(M)
    A <- K/mean(diag(K))
  }
  return(A)
}


D.mat <- function(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
                  n.core=1,shrink=FALSE,return.imputed=FALSE, ploidy=2, return.Xd=FALSE, 
                  method=3){
  
  #X <- apply(gg,2,function(x){y <- x; y[which(is.na(x))] <- mean(x, na.rm=TRUE); return(y)}); gg2[1:5,1:5]
  
  
  # ty <- apply(X, 2, function(x){length(table(x))})
  # vv <- which(ty > 2)
  # if(length(vv)==0){}
  
  # check if markers have dominant calls (0's)
  ty <- apply(X, 2, function(x){ 
    if(length(which(x==0) > 0)){
      return(1)}else{ # 1's are markers with dominance calls
        return(0)} # 0's are markers with only 1's or/and -1's
  }) 
  if(length(which(ty == 1)) < 2){ # there should be at least 2 markers with dominance calls to calculate the dominance relationship matrix
    cat("No heterozygous markers detected in the data. You might be using inbred lines.\nIf so, divide the markers in heterotic groups and do the kronecker product \namong A.mat's of the 2 groups to obtain a dominance relationship matrix\n")
    stop(call. = FALSE)
  }
  
  X2 <- X#[,vv]# only good markers with heterozygote plants
  # now transform 0 to 1's
  if(ploidy == 2){
    # Aa = 1 and AA|aa = 0
    X3 <- 1 - abs(X2)#apply(X2,2,function(x){y <- x; y[which(x == 1 | x==-1)] <- 0; y[which(x == 0)] <- 1; return(y)})
    if(return.Xd){
      X6 <- X3
    }else{
      if(method==1){
        X6 <- A.mat(X3, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                    n.core=n.core,shrink=shrink,return.imputed=return.imputed)
      }else if (method==2){
        
        M <- scale(X3, center = TRUE, scale = FALSE)
        K <- tcrossprod(M)
        
        bAlleleFrequency <- colMeans(X2+1)/2 # using original -1,0,1 matrix
        varHW <- sum((2 * bAlleleFrequency * (1 - bAlleleFrequency))^2) 
        
        X6 <- K/varHW
      }else if(method==3){
        #print("using")
        #X3 <- 1 - abs(X2)
        n <- dim(X2)[1]
        p <- colSums(X2+1)/(2*n) # from marker marix in 0,1,2 format
        q <- 1-p
        varHW <- sum(2*p*q * (1-(2*p*q)) )
        
        X3pq <- apply(X3, 1, function(x){ x - (2 * p * q)})
        X6 <- crossprod(X3pq)/varHW
        
      }
      
    }
  }else{
    X3 <- X2 - (ploidy/2)
    possible <- (-(ploidy/2):(ploidy/2))
    homo <- c(possible[1],possible[length(possible)])
    hete <- setdiff(possible, homo)
    X4 <- apply(X3,2,function(x){y <- x; y[which(x %in% hete)] <- 1; y[which(x %in% homo)] <- 0; return(y)})
    ty2 <- apply(X4, 2, function(x){length(table(x))})
    vv2 <- which(ty2 > 1)
    X5 <- X4[,vv2]# only good
    if(return.Xd){
      X6 <- X5
    }else{
      if(method==1){
        X6 <- A.mat(X5, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                    n.core=n.core,shrink=shrink,return.imputed=return.imputed)
      }else if (method==2){
        
        M <- scale(X5, center = TRUE, scale = FALSE)
        K <- tcrossprod(M)
        
        bAlleleFrequency <- colMeans(X5+1)/2 # using original -1,0,1 matrix
        varHW <- sum((2 * bAlleleFrequency * (1 - bAlleleFrequency))^2) 
        
        X6 <- K/varHW
      }else if(method==3){
        #X3 <- 1 - abs(X2)
        n <- dim(X5)[1]
        p <- colSums(X5+1)/(2*n) # from marker marix in 0,1,2 format
        q <- 1-p
        varHW <- sum(2*p*q * (1-(2*p*q)) )
        
        X5pq <- apply(X5, 1, function(x){ x - (2 * p * q)})
        X6 <- crossprod(X5pq)/varHW
        
      }
      
    }
    
  }
  return(X6)
}

E.mat <- function(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
                  n.core=1,shrink=FALSE,return.imputed=FALSE, type="A#A", ploidy=2){
  
  if(type == "A#A"){
    X4 <- A.mat(X, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                n.core=n.core,shrink=shrink,return.imputed=return.imputed, ploidy=ploidy)
    X5 <- hadamard.prod(X4,X4)
  }
  if(type == "A#D"){
    X4 <- A.mat(X, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                n.core=n.core,shrink=shrink,return.imputed=return.imputed, ploidy=ploidy)
    X4D <- D.mat(X, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                 n.core=n.core,shrink=shrink,return.imputed=return.imputed,ploidy=ploidy)
    X5 <- hadamard.prod(X4,X4D)
  }
  if(type == "D#D"){
    X4D <- D.mat(X, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                 n.core=n.core,shrink=shrink,return.imputed=return.imputed, ploidy=ploidy)
    X5 <- hadamard.prod(X4D,X4D)
  }
  return(X5)
}

pedtoK <- function(x, row="Row",column="Column",value="Ainverse", returnInverse=TRUE){
  
  # if(type=="asreml"){
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
    
    Ks <- as(K, Class = "sparseMatrix")
    if(returnInverse){
      Ksi <- solve(Ks)
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

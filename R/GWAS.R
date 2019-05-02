GWAS <- function(fixed, random, rcov, data, weights, 
                  iters=20, tolpar = 1e-03, tolparinv = 1e-06, 
                  init=NULL, constraints=NULL, method="NR", 
                  getPEV=TRUE,
                  na.method.X="exclude",
                  na.method.Y="exclude",
                  return.param=FALSE, 
                  date.warning=TRUE,
                  verbose=TRUE,
                  M=NULL, gTerm=NULL, n.PC = 0, min.MAF = 0.05, 
                  n.core=1, P3D = TRUE){
  
  if(is.null(gTerm)){
    stop("Please provide the name of the genetic term in the model in the gTerm argument", call. = FALSE)
  }
  
  qvalue <- function(p) {
    smooth.df = 3
    if (min(p) < 0 || max(p) > 1) {
      print("ERROR: p-values not in valid range.")
      return(0)
    }
    lambda = seq(0, 0.9, 0.05)
    m <- length(p)
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
    pi0 <- predict(spi0, x = max(lambda))$y
    pi0 <- min(pi0, 1)
    if (pi0 <= 0) {
      print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
      return(0)
    }
    u <- order(p)
    qvalue.rank <- function(x) {
      idx <- sort.list(x)
      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl
      return(tbl)
    }
    v <- qvalue.rank(p)
    qvalue <- pi0 * m * p/v
    qvalue[u[m]] <- min(qvalue[u[m]], 1)
    for (i in (m - 1):1) {
      qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                          1)
    }
    return(qvalue)
  }
  qq <- function(scores) {
    remove <- which(scores == 0)
    if (length(remove) > 0) {
      x <- sort(scores[-remove], decreasing = TRUE)
    }
    else {
      x <- sort(scores, decreasing = TRUE)
    }
    n <- length(x)
    unif.p <- -log10(ppoints(n))
    plot(unif.p, x, pch = 16, xlab = "Expected -log(p)", 
         ylab = "Observed -log(p)")
    lines(c(0, max(unif.p)), c(0, max(unif.p)), lty = 2)
  }
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }
  
  
  ## return all parameters for a mixed model
  res <- mmer(fixed=fixed, random=random, rcov=rcov, data=data, weights=weights, 
              iters=iters, tolpar=tolpar, tolparinv=tolparinv, 
              init=init, constraints=constraints, method=method, 
              getPEV=getPEV,
              na.method.X=na.method.X,
              na.method.Y=na.method.Y,
              return.param=TRUE, 
              date.warning=date.warning,
              verbose=verbose)
  # print(str(res))
  
  if(return.param){
    lastmodel <- res
  }else{
    ## fit initial models if P3D
    if(P3D){
      lastmodel <- .Call("_sommer_MNR",PACKAGE = "sommer",res[[1]], res[[2]],res[[3]],
                         res[[4]],res[[5]],res[[6]],res[[7]],res[[8]],res[[9]],
                         res[[10]],res[[11]],res[[12]], res[[13]],res[[14]],res[[15]],
                         TRUE)
      H2inv <- lastmodel$Vi
    }
    ## get names of random effects
    re_names <- unlist(res[[17]])
    # print(re_names)
    re_names <- gsub("\\(Intercept):","",re_names)
    
    gTermi <- which(re_names %in% gTerm)
    if(length(gTermi) < 1){
      stop("No match between the Gterm and the random effects present.")
    }
    
    nt <- ncol(as.matrix(res[[1]]))
    ntnames <- colnames(res[[1]])
    
    Ym <- matrix(t(scale(res[[1]])),ncol=1, byrow = FALSE) ##alternated
    
    Zb <- do.call(cbind,res[[4]][gTermi])
    Zm <- kronecker(Zb,diag(nt))
    
    if (n.PC > 0) {
      Kb <- do.call(adiag1,res[[5]][gTermi])
      eig.vec <- eigen(Kb)$vectors
      
      for(o in 1:length(res[[3]])){
        if(o == 1){
          Xm <- kronecker(res[[2]][[o]],res[[3]][[o]])
        }else{
          Xm <- cbind(Xm,kronecker(res[[2]][[o]],res[[3]][[o]])) 
        }
      }
      
      zbeig <- kronecker(Zb %*% eig.vec[,1:n.PC],diag(nt))
      Xm <- make.full(cbind(Xm,zbeig))
    }else {
      for(o in 1:length(res[[3]])){
        if(o == 1){
          Xm <- kronecker(res[[2]][[o]],res[[3]][[o]])
        }else{
          Xm <- cbind(Xm,kronecker(res[[2]][[o]],res[[3]][[o]])) 
        }
      }
      Xm <- make.full(Xm)
    }
    
    if(nrow(M) != ncol(Zb)){
      stop(paste("Marker matrix M needs to have same numbers of rows(",nrow(M),") than columns of the gTerm incidence matrix(",ncol(Zb),")."),call. = FALSE)
    }
    
    m <- ncol(M)
    colnamesM <- colnames(M)
    scorecalc <- function(i){
      
      Mi <- as.matrix(M[, i]);# colnames(Mi) <- colnamesM[i]
      Mi <- kronecker(Mi,diag(nt))
      freq <- mean(Mi + 1, na.rm = TRUE)/2
      MAF <- min(freq, 1 - freq)
      if (MAF >= min.MAF) {
        n2 <- length(Ym)
        y2 <- Ym
        Z2 <- Zm
        X3 <- cbind(Xm, Z2%*%Mi) 
        p <- ncol(X3)
        v1 <- 1
        v2 <- n2 - p
        if (!P3D) {## estimate varcomp for each marker 
          lastmodel <- .Call("_sommer_MNR",PACKAGE = "sommer",res[[1]], res[[2]],res[[3]],
                             res[[4]],res[[5]],res[[6]],res[[7]],res[[8]],res[[9]],
                             res[[10]],res[[11]],res[[12]], res[[13]],res[[14]],FALSE,
                             TRUE)
          H2inv <- lastmodel$Vi
        }
        W <- crossprod(X3, H2inv %*% X3)
        Winv <- try(solve(W), silent = TRUE)
        if (class(Winv) != "try-error") {
          xvy <- crossprod(X3, H2inv %*% y2)
          beta <- Winv %*% xvy # (XV-X)- XV-y
          # print(beta)
          resid <- y2 - X3 %*% beta # Y - XB
          s2 <- as.double(crossprod(resid, H2inv %*%resid))/v2 # eV-e/(n-p) = variance
          CovBeta <- s2 * Winv # eVe * B
          
          ps <- (ncol(Xm)+1):length(beta)
          res0 <- as.vector(beta)[ps]
          Fstat <- beta[ps]^2/diag(as.matrix(CovBeta[ps, ps]))
          x <- as.vector(v2/(v2 + v1 * Fstat))
          res1 <- -log10(pbeta(x, v2/2, v1/2))
          res2 <- Fstat
          
          SST <- t(y2)%*%y2
          SSR <- as.vector(t(beta)%*%xvy)
          SSM <- sum(y2)*mean(y2)
          # R-squared
          R2 = (SSR - SSM)/(SST - SSM)
          N <- length(y2); r<- ncol(X3)
          R2S = 1 - ( (N-1)*(1-R2)/(N-r) )
          res3 <- c(R2,R2S)
          
          myres <- matrix(c(res0,res1,res2,res3),ncol=1)
          colnames(myres) <- colnames(Mi)
          return(myres)
        }
      }
    }
    
    
    cat("Performing GWAS evaluation\n")
    if ((n.core > 1) & requireNamespace("parallel",quietly=TRUE)) {
      scores <- parallel::mclapply(as.list(1:m), function(x) {
        scorecalc(x)
      }, mc.cores = n.core)
    } else {
      scores <- lapply(as.list(1:m), function(x) {
        scorecalc(x)
      })
    }
    
    tokeep <- which(!unlist(lapply(scores, is.null)))
    
    scores <- do.call(cbind,scores[tokeep])
    
    # scores <- matrix(scores,nrow=(nt*3)+2,byrow = FALSE)
    
    mm <- as.data.frame(expand.grid(ntnames,c("beta","score","Fstat")))
    mm$tt <- paste(mm[,1],mm[,2])
    colnames(scores) <- colnames(M)[tokeep]
    # lastmodel$jkl <- colnames(M)[tokeep]
    rownames(scores) <- c(mm$tt,"R2","R2s")
    lastmodel$scores <- scores
    lastmodel$method <- method
    lastmodel$constraints <- res[[8]]
    # lastmodel$constraintsF <- res[[9]]
    class(lastmodel)<-c("mmergwas")
  }
  
  return(lastmodel)
}

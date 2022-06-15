# R implementation of EM, mmer and mmec have a c++ implementation
EM <- function(y,X=NULL,ZETA=NULL,R=NULL,iters=30,draw=TRUE,silent=FALSE, constraint=TRUE, init=NULL, forced=NULL, tolpar = 1e-04, tolparinv = 1e-06){
  convergence <- FALSE
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### fix that some initially were not specified
  if(is.null(X)){
    X <- matrix(1,nrow=length(y))
  }
  if(is.null(R)){
    R <- list(units=diag(length(y)))
  }
  y.or <- y
  X.or <- X
  ZETA.or <- ZETA
  R.or <- R
  if(is.null(names(ZETA))){
    varosss <- c(paste("u.",1:length(ZETA), sep=""))
  }else{
    varosss <- c(names(ZETA))
  }; varosssZ <- varosss
  if(is.null(names(R))){
    varosss <- c(varosss,paste("Res.",1:length(R),sep=""))
  }else{
    varosss <- c(varosss,names(R))
  }
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  #### %%%%%%%%%%%%%%%%%%%%%%%%
  
  #### now reduce the original inputs
  #y <- scale(y)
  good <- which(!is.na(y))
  y <- y[good]
  X <- as.matrix(X[good,])
  ZETA <- lapply(ZETA, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
  R <- lapply(R, function(x,good){x <- x[good,good]; return(x)}, good=good)
  #### get sizes of reduced
  qr <- qr(X)
  ranx <- length(y)-qr$rank # length(good)
  nz <- length(ZETA)
  nr <- length(R)
  nx <- dim(X)[2]
  dimzs <- lapply(ZETA,function(x){dim(as.matrix(x$Z))[2]})
  dimrs <- lapply(R,function(x){dim(x)[1]})
  N <- length(y)
  
  #### get the indexes for all effects
  sta <- 1
  ind1 <- numeric()
  for(u in 1:length(dimzs)){
    sta <- dimzs[[u]]+sta
    ind1[u] <- sta
  }; ind1<- c(1,ind1[-c(length(ind1))]);
  ind2 <- (c(ind1-1, sum(unlist(dimzs))))[-1]
  
  ind1g <- ind1 + nx ## indexes for Gi
  ind2g <- ind2 + nx ## indexes for Gi
  ind1 <- c(1,ind1 + nx) ## indexes for all including x
  ind2 <- c(nx,ind2 + nx) ## indexes for all including x
  
  #### initialize var components
  if(is.null(init)){
    varcom <- rep(var(y,na.rm=TRUE)/(nz+nr),nz+nr)
    names(varcom) <- c(rep("z",nz),rep("r",nr))
    varcomz <- varcom[which(names(varcom)=="z")]
    varcomr <- varcom[which(names(varcom)=="r")]
  }else{
    varcom <- init 
    if(length(init)!= (nr+nz)){
      stop("Please provide initial values for all variance components",call. = FALSE)
    }else{
      names(varcom) <- c(rep("z",nz),rep("r",nr))
      varcomz <- varcom[which(names(varcom)=="z")]
      varcomr <- varcom[which(names(varcom)=="r")]
    }
  }
  
  ll2=-10000000 # initial log-likelihood
  ll.stored <- ll2
  conv=0 # convergence
  wi=0 # counter
  taper <- rep(1, iters) # weighting parameter for updates
  #taper[1:3] <- c(.5,.7,1)#c(0.5, 0.7) # c(0.5, 0.7)
  
  #### get inverse of covariance and residual matrices
  Riw <-  lapply(R,function(x){solve(x)})
  Ki <- lapply(ZETA,function(x){
    findrank <- qr(x$K)$rank
    if(findrank < dim(x$K)[1]){
      return(solve(x$K + diag(tolparinv,dim(x$K)[1])))
    }else{
      return(solve(x$K))
    }
  })
  
  varso <- as.matrix(varcom)
  #### start
  if(!silent){
    count <- 0
    tot <- iters+1
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  if(!is.null(forced)){
    varcom <- forced
    if(length(init)!= (nr+nz)){
      stop("Please provide initial values for all variance components",call. = FALSE)
    }else{
      names(varcom) <- c(rep("z",nz),rep("r",nr))
      varcomz <- varcom[which(names(varcom)=="z")]
      varcomr <- varcom[which(names(varcom)=="r")]
    }
    iters=1
  }
  ################### =================================================================
  ################### =================================================================
  ################### =================================================================
  while (conv==0) { # ==================== START EM ALGORITHM =========================
    wi=wi+1
    if(!silent){
      count <- count + 1
    }
    ##############################
    ### for R * se as a direct sum
    Rse <- R[[1]]*0
    for(u in 1:length(R)){
      Rse <- Rse + R[[u]]
    }
    ##############################
    ## R inverse
    Rsei <- solve(Rse)#solve(Rse)
    ##############################
    ## G inverse
    varoz <- as.list(varcomz)  # variance components as list, no error included
    Gi <- lapply(as.list(c(1:nz)),function(x,Ko,v){
      oo=Ko[[x]]*as.numeric(varcomr[1]/(v[[x]])); return(oo)
    }, Ko=Ki, v=varoz) ## K*v(u)
    
    Zbind <- do.call(cbind,lapply(ZETA, function(x){x$Z}))
    XZbind <- as(cbind(X,Zbind), Class="sparseMatrix")
    CM <- t(XZbind) %*% Rsei %*% XZbind
    
    ## add Gi's to the C matrix
    for(u in 1:length(Gi)){
      rox <- ind1g[u]; cox <- ind2g[u]
      CM[rox:cox,rox:cox] <- CM[rox:cox,rox:cox] + Gi[[u]]
    }
    
    ## do gaussian eliminations (absorption)
    RHS <- t(XZbind) %*% Rsei %*% y
    ge <- solve(CM,RHS)
    ## invert the coefficient matrix
    CMi <- solve(CM)
    
    ##%%%%%%%%%%%%%%
    ### likelihood
    ##%%%%%%%%%%%%%%
    na <- sum(unlist(dimzs));na
    se <- varcomr
    #### Calculate ypy by building M
    left <- t(XZbind) %*% Rsei %*% y
    top <- t(left)
    corner <- t(y)%*%Rsei%*%y
    
    M <- rbind(cbind(corner,top),cbind(left,CM))
    M[1:4,1:4]
    vb <- try(chol(as(M, Class="sparseMatrix"),pivot=FALSE), silent=TRUE)
    if(is(vb, "try-error") ){ # class(vb)=="try-error"
      vb <- try(chol(as(M+diag(tolparinv,dim(M)[1]), Class="sparseMatrix")), silent=TRUE)
      
      if(is(vb, "try-error")){ # class(vb)=="try-error"
        ypy <- 1#(L[1,1]^2); ypy
        logdc <- 2#2* sum(log(diag(L)[-1])); logdc
      }else{
        L <- t(vb); L[1:4,1:4]
        ypy <- (L[1,1]^2); ypy
        logdc <- 2* sum(log(diag(L)[-1])); logdc
      }
    }else{
      L <- t(vb); L[1:4,1:4]
      ypy <- (L[1,1]^2); ypy
      logdc <- 2* sum(log(diag(L)[-1])); logdc
    }
    
    
    logda <- sum(unlist(lapply(ZETA,function(x){determinant(x$K, logarithm = TRUE)$modulus[[1]]})))
    nlogsu <- log(varcomz)*unlist(dimzs) # n*log(su)
    
    
    ll <- - 0.5 * ( ((N - ranx - na)*log(se)) + logdc + logda + sum(nlogsu) +  ypy);ll
    if (abs(ll-ll2) < tolpar | wi == iters ){ ## CONVERGENCE, not sure if should be absolute value or not
      conv=1
      if (abs(ll-ll2) < tolpar){
        convergence <- TRUE
      }
    }
    ll2=(ll) # replace the initial logLik value for the value reached
    ll.stored <- c(ll.stored,ll2)
    ##%%%%%%%%%%%%%%
    ##%%%%%%%%%%%%%%
    
    if(!silent){
      setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
    }
    
    ##%%%%%%%%%%%%%%
    ### UPDATE PARAMETERS
    ##%%%%%%%%%%%%%%
    ## extract BLUPs and BLUEs
    ulist <-list()
    for(k in 1:length(ind1g)){
      ulist[[k]] <- as.matrix(ge@x[ind1g[k]:ind2g[k]])
    }
    b <- as.matrix(ge@x[1:nx])
    Xb <- (X%*%b)
    
    Zulist <- list()
    Welist <- list()
    for(k in 1:length(ind1g)){
      Zulist[[k]] <- ZETA[[k]]$Z %*% ulist[[k]]
      Welist[[k]] <- Zulist[[k]]/varcomz[k]
    }
    zu <- do.call(cbind,Zulist); zu <- rowSums(zu)
    ##############################
    ## new estimate for error variance
    ## y'y - new.y'y / (N - nx)
    ##############################
    now <- 0 
    
    for(f in 1:length(ind1g)){ # y'Xb + y'Zu
      now <- now + crossprod(as.matrix(Zulist[[f]]),y)
    }
    b <- as.matrix(ge@x[ind1[1]:ind2[1]])
    now <- now + crossprod(X%*%b,y)
    
    varcomr[1] <- ( (t(y)%*%as.matrix(Rsei)%*%y) - now ) / (length(y)-nx)
    
    ##############################
    ## new estimates for variance components
    ##############################
    for(k in 1:length(ind1g)){ # adjust var comps except ERROR VARIANCE
      Kinv <- Ki[[k]]
      nk <- nrow(Kinv)
      rox <- ind1g[k]; cox <- ind2g[k]
      varcomz[k] <- ((t(ulist[[k]]) %*% Kinv  %*% ulist[[k]] ) + (as.numeric(varcomr[1])*sum(diag(Kinv%*%CMi[rox:cox,rox:cox]))))/nk 
    }
    
    #print(c(varcomz,varcomr))
    varcom <- c(varcomz,varcomr)
    zeros <- which(varcom <= 0)
    
    if(length(zeros)>0 & constraint){
      varcom[zeros] <- varcomr[1] * (1.011929e-07^2)
    }
    ## move back to the z and r separation
    names(varcom) <- c(rep("z",nz),rep("r",nr))
    varcomz <- varcom[which(names(varcom)=="z")]
    varcomr <- varcom[which(names(varcom)=="r")]
    #y <- Xb+zu
    varso <- cbind(varso,as.matrix(varcom))# just to keep track
  } ################# ======================== END OF ALGORITHM =======================
  ################### =================================================================
  ################### =================================================================
  ################### =================================================================
  e <- y - (Xb+zu)
  
  
  We <- (e)/varcomr[1]
  B <- cbind(do.call(cbind,Welist),We)
  left <- t(XZbind) %*% Rsei %*% B #left of the M matrix
  top <- t(left) # top of M
  corner <- t(B)%*%Rsei%*%B # left corner of M
  M <- rbind(cbind(corner,top),cbind(left,CM))
  vb <- try(chol(as(M, Class="sparseMatrix")), silent = TRUE)
  if( is(vb, "try-error") ){ # class(vb)=="try-error"
    aii <- diag(nz+nr)
  }else{
    L <- t(vb); #L[1:4,1:4]
    ai <- L[1:(nz+nr),1:(nz+nr)]^2;ai
    ai <- as.matrix(ai)#*100#(var(y,na.rm = TRUE)) #originally should not be multiplied
    ai[upper.tri(ai)] <- t(ai)[upper.tri(ai)]
    aii <- solve(ai)
  }
  
  
  ## plot likelihood
  if(draw){
    layout(matrix(1:2,2,1))
    plot(y=ll.stored[-1], x=1:length(ll.stored[-1]), type="l", main="Log-likelihood", xlab="iteration",ylab="Log-likelihood")
    ## plot varcomp
    for(l in 1:dim(varso)[1]){
      #lb <- min(unlist(varso),na.rm = TRUE)
      ub <- max(unlist(varso),na.rm = TRUE)
      if(l==1){plot(varso[1,],ylim=c(0,ub),type="l", main="MME-EM results",ylab="varcomp", xlab="iteration")}else{
        lines(varso[l,],x=1:dim(varso)[2], col=l)
      }
    }
  }
  ## ======================== ##
  ## ======================== ##
  ## PROVIDE EXTRA PARAMETERS
  ## ======================== ##
  ## ======================== ##
  
  ## variance components
  out1 <- as.matrix(varcom); rownames(out1) <- varosss; colnames(out1) <- "component"
  ## inverse of the phenotypic variance (ZKZ'+R)-
  Vinv <- XZbind %*% CMi %*% t(XZbind)
  ## BLUPs
  names(ulist) <- varosssZ
  for(f in 1:length(ZETA)){
    rownames(ulist[[f]]) <- colnames(ZETA[[f]]$Z)
  }
  ## VarBLUPs
  Var.u <- list()
  for(f in 1:length(ind1g)){
    Var.u[[f]] <- diag(CMi[ind1g[f]:ind2g[f],ind1g[f]:ind2g[f]])
    names(Var.u[[f]]) <- colnames(ZETA[[f]]$Z)
  }
  ## betas
  rownames(b) <- colnames(X)
  ## var.betas
  xvxi <- CMi[ind1[1]:ind2[1],ind1[1]:ind2[1]]
  rownames(xvxi) <- colnames(xvxi) <- colnames(X)
  ## cond. residuals
  #e
  ## residuals
  ee <- y - (Xb)
  ## log likelihood
  ll <- as.vector(ll)
  ## AIC and BIC
  AIC = as.vector((-2 * ll ) + ( 2 * nx))
  BIC = as.vector((-2 * ll ) + ( log(length(y)) * nx))
  ## fish.inv
  # aii
  ## Ksp
  Ksp <- do.call(adiag1,lapply(ZETA,function(x){x$K}))
  ## fitted values only for non-missing data
  fit0 <- Xb + zu
  
  ## ======================== ##
  ## ======================== ##
  ## 2. PROVIDE EXTRA PARAMETERS
  ## using original data
  ## ======================== ##
  ## ======================== ##
  Xor.b <- X.or%*%b 
  Zor.u <-  lapply(as.list(1:nz),function(x,z,u){z[[x]]$Z %*% ulist[[x]]},z=ZETA.or,u=ulist)
  Zor.u2 <- do.call(cbind,Zor.u); Zor.u2 <- rowSums(Zor.u2)
  fit1 <- Xor.b + Zor.u2 # fitted values
  
  if(!silent){
    setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
  }
  
  ### let user knows the constraints
  out1 <- as.data.frame(out1)
  out1[,"constraint"] <- "Positive"
  if(length(zeros)>0){
    out1[zeros,"constraint"] <- "Boundary"
    out1[zeros,1] <- 0
  }
  
  sigma.scaled <- out1[,1]/var(y,na.rm=TRUE)
  
  res <- list(var.comp=out1, V.inv = Vinv, u.hat=ulist, Var.u.hat=Var.u, 
              #PEV.u.hat=PEV.u, 
              beta.hat=b, Var.beta.hat=xvxi, residuals=ee, cond.residuals=e,
              LL=ll, sigma.scaled=sigma.scaled,
              AIC=AIC, BIC=BIC, fish.inv=aii,fitted.y.good=fit0, 
              X=X.or, Z=Zbind, K=Ksp, ZETA=ZETA,
              ## go back to original data
              fitted.y=fit1, fitted.u=Zor.u2, 
              forced=forced, convergence=convergence,
              CM=CM, CMi=CMi)
  
  layout(matrix(1,1,1))
  return(res)
}

# R implementation of EMMA, mmer and mmec have a c++ implementation
MEMMA <- function (Y, X=NULL, ZETA=NULL, tolpar = 1e-06, tolparinv = 1e-06, check.model=TRUE, silent=TRUE) {
  
  Y <- as.matrix(Y)
  if(is.null(colnames(Y))){
    colnames(Y) <- paste("T",1:ncol(Y),sep=".")
  }
  if(is.null(X)){
    X <- matrix(1,nrow=dim(Y)[1])
  }
  respo <- colnames(Y)
  
  if(check.model){ # if needs to be checked, else just skip
    ZETA <- lapply(ZETA, function(x){
      if(length(x) == 1){
        provided <- names(x)
        if(provided == "Z"){
          y <- list(Z=x[[1]],K=diag(dim(x[[1]])[2]))
        }else if(provided == "K"){
          y <- list(Z=diag(length(y)),K=x[[1]])
        }else{
          stop(call.=FALSE)
          cat("Names of matrices provided can only be 'Z' or 'K', the names you provided don't match the arguments required")
          jkl <- c(23,18,9,20,20,5,14, NA,2,25,NA,7,9,15,22,1,14,14,25,NA,3,15,22,1,18,18,21,2,9,1,19)
          oh.yeah <- paste(letters[jkl],collapse = "")
        }
      }else{y <- x}; 
      return(y)
    })
  }
  
  havetobe <- apply(Y,2,is.numeric)
  if(length(which(havetobe)) != dim(Y)[2]){
    stop("The response variables need to be numeric\n", call.=FALSE)
  }
  Y <-apply(Y,2, function(x){vv<-which(is.na(x)); if(length(vv)>0){x[vv]<-mean(x,na.rm=TRUE)};return(x)})
  
  Zlist <- lapply(ZETA, function(x){x$Z})
  Klist <- lapply(ZETA, function(x){x$K})
  Zs <- do.call("cbind", Zlist)
  Ks <- do.call("adiag1", Klist)
  
  K <- Zs%*%Ks%*%t(Zs)
  Z <- diag(dim(K)[1])
  
  X <- t(X)
  Y <- t(Y)
  
  dim(Z);dim(X);dim(Y); dim(K)
  #Z <- t(Z)
  ECM1 <- function(ytl, xtl, Vgt, Vet, Bt, deltal) {
    ## deltal(de) is the eigen decomposition of ZKZ'
    ##  V = de * Vg' + Ve' 
    Vlt = deltal * Vgt + Vet
    ## Vinv (add some noise to make sure is invertible)
    invVlt <- solve(Vlt + tolparinv * diag(d))
    ## Vlt = V
    ## gtl = Vg'*Vinv*(y-Xb)
    ## Sigmalt = de*Vg' - de*(Vg'*Vinv*(de*Vg'))
    return(list(Vlt = Vlt, gtl = deltal * Vgt %*% invVlt %*% 
                  (ytl - Bt %*% xtl), Sigmalt = deltal * Vgt - deltal * 
                  Vgt %*% invVlt %*% (deltal * Vgt)))
  }
  # for each trait extract response and do eigen decomposition of eigen
  wrapperECM1 <- function(l) { 
    ytl <- Yt[, l]
    xtl <- Xt[, l]
    deltal <- eigZKZt$values[l]
    return(ECM1(ytl = ytl, xtl = xtl, Vgt = Vgt, Vet = Vet, 
                Bt = Bt, deltal = deltal))
  }
  # genetic variance
  Vgfunc <- function(l) {
    # gtl = Vg'*Vinv*(y-Xb), then gtl*gtl
    Vgl <- tcrossprod(outfromECM1[[l]]$gtl)
    # (1/n) * (1/eigvalues) * [[[ Vg'*Vinv*(y-Xb).*.Vg'*Vinv*(y-Xb) ]]] * [[[ de*Vg' - de*(Vg'*Vinv*(de*Vg')) ]]]
    return((1/n) * (1/eigZKZt$values[l]) * (Vgl + outfromECM1[[l]]$Sigmalt))
  }
  Vefunc <- function(l) {
    ## error' trait l
    ## Y' - B'X' -  Vg'*Vinv*(y-Xb)
    etl <- Yt[, l] - Bt %*% Xt[, l] - outfromECM1[[l]]$gtl
    ## return (1/n) * e'e + [[de*Vg' - de*(Vg'*Vinv*(de*Vg'))]]
    return((1/n) * ((tcrossprod(etl) + outfromECM1[[l]]$Sigmalt)))
  }
  if (sum(is.na(Y)) == 0) {
    N <- nrow(K)
    KZt <- tcrossprod(K, Z)
    ZKZt <- Z %*% KZt
    eigZKZt = eigen(ZKZt)
    n <- nrow(ZKZt)
    d <- nrow(Y)
    Yt = Y %*% eigZKZt$vectors
    Xt = X %*% eigZKZt$vectors
    Vgt = cov(t(Y))/2
    Vet = cov(t(Y))/2
    XttinvXtXtt <- t(Xt) %*% solve(tcrossprod(Xt))
    Bt <- Yt %*% XttinvXtXtt
    Vetm1 <- Vet
    repeat {
      outfromECM1 <- lapply(1:n, wrapperECM1)
      Vetm1 <- Vet
      Gt = sapply(outfromECM1, function(x) {
        cbind(x$gtl)
      })
      Bt = (Yt - Gt) %*% XttinvXtXtt
      listVgts <- lapply(1:n, Vgfunc)
      Vgt <- Reduce("+", listVgts)
      listVets <- lapply(1:n, Vefunc)
      Vet <- Reduce("+", listVets)
      convnum <- abs(sum(diag(Vet - Vetm1)))/abs(sum(diag(Vetm1)))
      convcond <- tryCatch({
        convnum < tolpar
      }, error = function(e) {
        return(FALSE)
      })
      if (convcond) {
        break
      }
    }
    ## V inverse
    HobsInv <- solve(kronecker(ZKZt, Vgt) + kronecker(diag(n),Vet) + tolparinv * diag(d * n))
    
    #print(dim(Y));print(dim(Bt));print(dim(X))
    ehat <- matrix(Y - Bt %*% X, ncol = 1, byrow = F) # residuals
    HobsInve <- HobsInv %*% ehat # V- (Y-XB)
    varvecG <- kronecker(K, Vgt) # G
    ## u.hat GZ'V-(Y-XB)
    gpred <- varvecG %*% (kronecker(t(Z), diag(d))) %*% HobsInve
    Gpred <- matrix(gpred, nrow = nrow(Y), byrow = F) # u.hat as matrix
    colnames(Gpred) <- rownames(K)
    Xforvec <- (kronecker(t(X), diag(d)))
    Zforvec <- (kronecker((Z), diag(d)))
    ZKforvec <- Zforvec %*% varvecG
    
    xvx <- crossprod(Xforvec,HobsInv %*% Xforvec)
    P <- HobsInv - HobsInv %*% Xforvec %*% solve(xvx, crossprod(Xforvec, HobsInv))
    ddv <- determinant(HobsInv, logarithm = TRUE)$modulus[[1]]
    Yvect <- as.matrix(as.vector(as.matrix(Y))) #dim(Y.or2)
    
    ytPy <- t(Yvect)%*%(P%*%(Yvect))
    llik=as.numeric(-0.5*((ddv)+determinant(solve(xvx), logarithm = TRUE)$modulus[[1]]+ytPy)) # log likelihood, problem
    
    varGhat <- crossprod(ZKforvec, P) %*% ZKforvec
    if (!exists("P")) {
      P <- HobsInv - HobsInv %*% Xforvec %*% solve(crossprod(Xforvec, 
                                                             HobsInv %*% Xforvec), crossprod(Xforvec, HobsInv))
    }
    PEVGhat <- varvecG - varGhat
    varBhat <- solve(crossprod(Xforvec, HobsInv %*% Xforvec))
    
    ######## AIC BIC
    AIC = as.vector((-2 * llik ) + ( 2 * dim(X)[1]))
    BIC = as.vector((-2 * llik ) + ( log(dim(as.matrix(Y))[2]) * dim(X)[1]))
    
    #print(varvecG)
    ehat <- t(Y) - t(X)%*%t(Bt) # residuals = Y - XB
    
    cond.ehat <- t(Y) - ( (t(X)%*%t(Bt)) + (Z %*% t(Gpred)) ) # cond.residuals = Y - (XB+Zu)
    
    fitted <- ( (t(X)%*%t(Bt)) + (Z %*% t(Gpred)) )
    
    fitted.u <-  Z %*% t(Gpred) 
    
    sigma <- list(Vu=Vgt, Ve=Vet)
    
    dimos <- lapply(ZETA, function(x){dim(x$Z)})
    
    
    u.hat <- t(Gpred)#unique(u.hat0)
    colnames(u.hat) <- respo
    #u.hat0 <- (t(Gpred))
    Z1 <- Zlist[[1]]
    namesZ1 <- colnames(Z1)
    if(!is.null(namesZ1)){
      rownames(u.hat) <- apply(Z1,1,function(x,y){paste(y[which(x==1)], collapse=".")},y=colnames(Z1))
    }
    
    #colnames(u.hat0) <- respo
    
    return(list(var.comp=sigma, V.inv=HobsInv, u.hat = u.hat , LL=llik, AIC=AIC,BIC=BIC,
                Var.u.hat = (varGhat), beta.hat = t(Bt),  Var.beta.hat = (varBhat), 
                PEV.u.hat = (PEVGhat), residuals=ehat, cond.residuals=cond.ehat,
                fitted.y=fitted, fitted.u=fitted.u, Z=Z, K=K, dimos=dimos, ZETA=ZETA,
                method="EMMAM", convergence=TRUE)) # XsqtestB = XsqtestB, pvalB = p.adjBhat, XsqtestG = XsqtestG,  pvalG = p.adjGhat,
  }
}

# R implementation of AI, mmer and mmec have a c++ implementation
# TO DO
# implement AR1 and FA models
## PARAMETER DETAILS
## X is the design matrix
## Z is a list of lists each element contains the Z matrices required for the covariance structure specified for a random effect
## Ai is a list with the inverses of the relationship matrix for each random effect
## y is the response variable
## S is the list of residual matrices
## H is the matrix of weights. This will be squared via the cholesky decomposition and apply to the residual matrices
## nIters number of REML iterations 
## tolParConvLL rule for stoping the optimization problem, difference in log-likelihood between the current and past iteration
## tolParConvNorm rule for stoping the optimization problem, difference in norms
## tolParInv value to add to the diagonals of a matrix that cannot be inverted because is not positive-definite
## theta is the initial values for the vc (matrices should be symmetric).
## thetaC is the constraints for vc: 1 positive, 2 unconstrained, 3 fixed
## thetaF is the dataframe indicating the fixed constraints as x times another vc, rows indicate the variance components, columns the scale parameters (other VC plus additional ones preferred)
## addScaleParam any additional scale parameter to be included when applying constraints in thetaF
## weightInfEMv is the vector to be put in a diagonal matrix (a list with as many matrices as iterations) representing the weight assigned to the EM information matrix
## weightInfMat is a vector of weights to the information matrix for the operation delta = I- * dLu/dLx # unstructured models may require less weight to the information matrix

AI <- function(X=NULL,Z=NULL, Zind=NULL, Ai=NULL,y=NULL,S=NULL, H=NULL,
                    nIters=80, tolParConvLL=1e-4, tolParConvNorm=.05, tolParInv=1e-6,
                    theta=NULL, thetaC=NULL, thetaF=NULL,addScaleParam=NULL,
                    weightInfEMv=NULL,weightInfMat=NULL){
  
  
  thetaCs <- lapply(thetaC, function(x){x[lower.tri(x)] <- t(x)[lower.tri(x)];return(x)})
  y <- as.matrix(y) # ensure y is a matrix IGNORE IN RCPP
  logDetA <- lapply(Ai, function(x){-determinant(x,logarithm = TRUE)$modulus}) # times minus one because we use the inverse
  thetaList <- list() # store thetas (variance components)
  llik <- numeric() # store log likellihood values
  
  # for each random effect identify how many levels exist and assign an identifier for #of effect
  last=ncol(X)
  partitions <- list()
  zsAva <- unique(Zind)
  
  if(!is.null(Z)){ # if random effects exist assign an identifier for the random effects
    for(j in 1:length(zsAva)){ # for each random effect
      use <- which(Zind == j)
      Nus = unlist(lapply(Z[use], function(x){ncol(x)})) # ncol(X) and ncol(y's)
      end <- Nus # dummy
      for(i in 1:length(Nus)){end[i]<- sum(Nus[1:i])} # real end (#of effect)
      start <- end-Nus+1 # start id for each random effect
      start=start+last; end=end+last
      partitions[[j]] <- cbind(start,end,j) # matrix telling us where each random effect starts and where it ends
      last=end[length(end)] # new last for the j'th random effect
    }; partitions
  }
  
  # total number of variance components per random effect (e.g., 1, 2, 3, ...)
  Nvc2 <- unlist(lapply(thetaC,function(x){ # for each thetaC matrix 
    s1 <- x[upper.tri(x, diag = TRUE)] # extract upper triangular and diagonal
    return(length(which(s1 > 0))) # how many values are different than zero
  }))
  
  # assign a start and an end to each covariance structure using the #of VC 
  endvc2 <- Nvc2
  for(i in 1:length(Nvc2)){endvc2[i]<- sum(Nvc2[1:i])}
  startvc2 <- endvc2-Nvc2+1
  
  # create a VECTOR with the contraints
  thetaCUnlisted= unlist(lapply(as.list(1:length(thetaC)), function(x){
    return(thetaC[[x]][which(thetaC[[x]]>0)])
  }))
  
  Nb = ncol(X) # number of fixed effects
  Nr = nrow(y) # number of records or observations
  Nu=Nus=0
  if(!is.null(Z)){ # only if random effects exist calculate the number of effects to estimate in the random part
    Nus = unlist(lapply(partitions, function(x){sum(x[1,2]-x[1,1]+1)})) # should we only consider the ncols from the first matrix in the case when covariance components exist?
    Nu = sum(Nus) # total number of random effects to be estimated
  }
  
  weightInfEMl <- list()
  for(s in 1:nIters){
    weightInfEMl[[s]] <- diag( sum(Nvc2) ) * weightInfEMv[s] # rep( list( diag( sum(Nvc2) ) * weightInfEMv[s] ), nIters) # a 0 weights for EM in general
  }
  
  # objects for storing
  finalThetaList <- list() # list to store the variance components for each iteration
  normChange1 <- normChange3 <- numeric() # vector to store the change 
  
  percChangeMat <- matrix(NA,nrow=length(thetaCUnlisted),ncol=nIters)
  changeNorm <- matrix(NA,3,nIters)
  toBoundary <- matrix(0,nIters,length(thetaCUnlisted))
  sumToBoundary <- numeric(length(thetaCUnlisted))
  ###########################
  # START ITERATIVE ALGORITHM
  ###########################
  for(iIter in 1:nIters){ # iIter=1
    print(paste("iteration",iIter,"-", strsplit(as.character(Sys.time())," ")[[1]][2])) # print iteration
    thetaList[[iIter]] <- do.call(adiag1,theta) # storo VC for the ith iteration
    
    ep = length(theta) # position of the error VC
    if(iIter==1){ # unlist theta (all VCs) to calculate at the end the % change of vc
      finalTheta <- unlistThetaWithThetaC(theta,thetaC)
    }
    ###########################
    # 1) absorption of m onto y
    # PAPER FORMULA from Jensen and Madsen 1997, Gilmour et al., 1995
    # expand coefficient matrix (C) to have the response variable
    # M = W' Ri W # with W = [X Z y] 
    #
    #     [X'RiX  X'RiZ     X'Riy ]
    # M = [Z'RiX  Z'RiZ+Gi  Z'Riy ]
    #     [y'RiX  y'RiZ     y'Riy ]
    #
    # where Gi = Ai*(s2e/s2u) = (A*s2u)*s2e = kronecker(Ai,solve(s2u))
    #
    # lambda = solve(theta) # inverse of var-covar matrices
    # MChol = chol(M)
    # yPy = MChol[n,n] # where n is the last element of the matrix
    # logDetC = 2 * E log(diag(MChol))
    ###########################
    if(!isDiagonal(H)){
      Hs <- chol(H)
    }
    Rij <- Rij.inv <- list()
    usedResidual <- which(thetaC[[length(thetaC)]] > 0)
    for(i in 1:length(S)){
      uri <- usedResidual[i]
      Rij[[i]] <- S[[i]]*(theta[[length(theta)]][uri]) 
      Rij.inv[[i]] <- S[[i]]*(1/theta[[length(theta)]][uri]) # do differently if inverting a complex S
    }  # R = Reduce("+",Rij)
    Ri = Reduce("+",Rij.inv)
    if(!isDiagonal(H)){
      Ri <- Hs %*% Ri %*% t(Hs)
    }
    
    M0 <- mmeFormation(X=X,Z=Z,Ri=Ri,y=y) # make a control to function if Ri is not provided
    M <- M0$M # all MMEs
    W <- M0$W # W = [X Z y] 
    
    lambda <- list() # store inverses of thetas (VCs matrices)
    GI <- list() # store inverses of covariance matrices for random effects
    
    if(!is.null(Z)){ # if random effects exist calculate inverses of VC matrices and Gi to add to MME
      for(iR in 1:length(zsAva)){ # for each random effect
        lambda[[iR]] = (solve(theta[[iR]]))
      }
      for(iR in 1:length(zsAva)){  # for each random effect
        # build Gi
        GI[[iR]] = kronecker(lambda[[iR]],Ai[[iR]]) # Ai ** lambda
        # add it to the MMEs 
        partitionsP=partitions[[iR]];
        iRpartition=partitionsP[1,1]:partitionsP[nrow(partitionsP),2]
        M[iRpartition,iRpartition] <- M[iRpartition,iRpartition] + GI[[iR]]
      }
    }
    
    MChol=try(chol(M,pivot=F), silent = TRUE) # Cholesky decomposition of C expanded by Y
    if(is(MChol, "try-error")){
      MChol=try(chol(M+diag(tolParInv,nrow(M)),pivot=F), silent = FALSE)
    }
    
    yPy = ((MChol[nrow(MChol),ncol(MChol)])^2) # last element of Cii = yPy
    logDetC = (2 * sum(log(diag(MChol)[-ncol(MChol)])) ) # sum(log(diag(Cii)))
    ###########################
    # 1.1) calculate the log-likelihood 
    # PAPER FORMULA (Lee and Van der Werf, 2006)    #
    # LL = -0.5 [((Nr-Nb-Nu-...)*ln(s2e)) - ln|C| + ln|Au| + ... + (Nu*ln(s2u)) + ... + y'Py ]
    # PAPER FORMULA (Jensen and Madsen, 1997)
    # LL = -0.5 [ln|C| + ln|R| + (ln|A.u| +ln|theta.u|) + ... + y'Py ]
    # where | | is the determinant of a matrix
    #       A.u is the pure relationship matrix for the uth random effect
    #       theta.u is the vc matrix for the uth random effect
    ###########################
    llikp <- 0
    if(!is.null(Z)){ # if random effects exist calculate ln|A.u| + ln|theta.u|
      for(iR in 1:length(zsAva)){ # for each random effect
        llikp = llikp + ((Nus[[iR]])*determinant(theta[[iR]],logarithm = TRUE)$modulus) + logDetA[[iR]]
      }
    }
    logDetR <- Nr * determinant(theta[[ep]],logarithm = TRUE)$modulus
    llik[iIter] = -.5 * ( llikp + logDetC + logDetR + yPy )
    
    ###########################
    # 2) backsubstitute to get b and u (CORRECT)
    # use the results from the absorption to obtain BLUE & BLUPs
    # b = backsolve(MChol[,rest],MChol[,last])
    ###########################
    bu=backsolve(MChol[,-ncol(MChol)],MChol[,ncol(MChol)])
    b=bu[1:Nb]
    u=NULL
    if(!is.null(Z)){ # if random effects exist extract u
      u=bu[(Nb+1):length(bu)]
    }
    ###########################
    # 3) calculate Wu (working variates)
    # PAPER FORMULA (Notes on Estimation of Genetic Parameters from Van der Werf)
    # wu = Zu/s2u; we = e/s2e
    # PAPER FORMULA (Jensen and Madsen, 1997)
    # U = [u1 | u2 | ... | ui]
    # US = U * lambda
    # Wu.ii = Zui*USi # for variance component
    # Wu.ij = Zui*USj + Zuj*USi # for covariance component
    # Wr.j = Rj * Rinv * e  # for residual variance component
    ###########################
    
    
    Zu <- Wu <- uu <- list()
    uAllList <- list()
    U.SinvList <- list()
    if(!is.null(Z)){ # if random effects exist
      for(iR in 1:length(zsAva)){ # for each random effect
        partitionsP = partitions[[iR]]
        uList <- list()
        for(irow in 1:nrow(partitionsP)){
          usedPartition = partitionsP[irow,1]:partitionsP[irow,2]
          uList[[irow]] = bu[usedPartition]#*4.762357 # u
        }
        U <- do.call(cbind,lapply(uList,as.matrix))
        uAllList[[iR]] <- U
        # [a || m] [s2a || sam] = [s2a a + sam m  || sam a + s2m m]
        #          [sam || s2m]
        U.Sinv <- U %*% lambda[[iR]]
        U.SinvList[[iR]] <- U.Sinv
        
        useZind <- which(Zind == iR)
        WuiR <- ZuiR <- list();counter=1
        # vcs2 <- numeric()
        for(irow in 1:nrow(lambda[[iR]])){
          for(icol in irow:ncol(lambda[[iR]])){
            # print(paste("irow",irow, "icol", icol))
            if(thetaC[[iR]][irow,icol] > 0){# if vc has to be estimated
              if(irow == icol){ # var comp
                ## Wu
                WuiR[[counter]] <- Z[[useZind[irow]]] %*% as.matrix(U.Sinv[,irow])
                ## Zu
                ZuiR[[counter]] <- Z[[useZind[irow]]] %*% as.matrix(U[,irow])
              }else{ # cov comp
                ## Wu
                WuiR[[counter]] <- ( Z[[useZind[icol]]] %*% U.Sinv[,irow] - Z[[useZind[irow]]] %*% U.Sinv[,icol] ) #* 2
                # WuiR[[counter]] <- ( Z[[useZind[icol]]] %*% U.Sinv[,irow] + Z[[useZind[irow]]] %*% U.Sinv[,icol] ) / 2
                ## Zu
                # ZuiR[[counter]] <- ( Z[[iR]][[icol]] %*% U[,irow] + Z[[iR]][[irow]] %*% U[,icol] )
              }
              counter=counter+1
            } # if vc has to be estimated
            
          }# for icol in 1row:nolc()
        }# for irow in 1:nrow()
        Wu[[iR]] = do.call(cbind, WuiR)
        Zu[[iR]] = do.call(cbind, ZuiR)
      }
    }
    
    Xb = X%*%b
    if(is.null(Z)){ 
      e = y - Xb # e = y - (Xb+Zu1+Zu2)
      iR=0
    }else{
      e = y - (Xb + Reduce("+",lapply(Zu,function(x){matrix(apply(as.matrix(x),1,sum))}))) # e = y - (Xb+Zu1+Zu2)
    }
    
    # useRes <- which(thetaC[[ep]] > 0)
    for(iS in 1:length(S)){ # for each residual effects
      Wu[[iR+iS]] = (S[[iS]]%*%Ri%*%e) # Wu.r = S * Rinv * e
    }
    ###########################
    # 4) absorption of m onto W (2 VAR, 1 COV)
    # we had to change the avInf to avInf/sigmas
    # PAPER FORMULA (Smith, 1995) Differentiation of the Cholesky Algorithm
    # avInf.ij = ((chol(M))[n,n])^2 # the square of the last diagonal element of the cholesky factorization
    # where:
    #        [X'RiX  X'RiZ     X'Riwj ]
    # M.Wu = [Z'RiX  Z'RiZ+Gi  Z'Riwj ]
    #        [wk'RiX wk'RiZ    wk'Riwj]
    # where:
    # wi: working variate i
    # wk: working variate k
    # Ri: is R inverse
    # and the the part corresponding to X and Z is the coefficient matrix C
    # AI = (M.Wu.chol)^2
    ###########################
    
    Wx = as.matrix(do.call(cbind,Wu))# [w1 w2 w3 ... wi]
    avInf = matrix(NA, nrow=ncol(Wx), ncol=ncol(Wx)) # AI
    
    C=M[-nrow(M),-ncol(M)] # remove y portion from M to get the coefficient matrix C # C11 upper left 
    XWj.ZWj = t(W)%*%Ri%*%Wx # [X'Riwj Z'Riwj]' # C12 upper right 
    WiX.WiZ = t(Wx)%*%Ri%*% W # [wk'RiX wk'RiZ] # C21 lower left
    WiWj= t(Wx)%*%Ri%*%Wx # wk'Riwj # C22 lower right
    MWu=rbind( cbind(C,(XWj.ZWj)) , cbind((WiX.WiZ),WiWj) ) # rbind( cbind(C11 C12), cbind(C21 C22) )
    
    MWuChol=try(chol(as.matrix(MWu),pivot=F), silent = TRUE) # Cholesky decomposition of C expanded by Wu
    if(is(MWuChol, "try-error")){
      MWuChol=try(chol(MWu+diag(tolParInv,nrow(MWu)),pivot=F), silent = FALSE)
    }
    
    nTheta = endvc2[length(endvc2)] # total number of variance components
    last=nrow(MWuChol); take = (last-nTheta+1):last
    avInf = MWuChol[take,take]
    avInf = t(avInf) %*% avInf
    # avInf = as.matrix(MWuChol[take,take]^2)
    # avInf[lower.tri(avInf)] <- t(avInf)[lower.tri(avInf)] # lower and upper are not equal (upper seems to be the right one to use)
    avInfInv <- try(solve(avInf), silent = TRUE)
    if(is(avInfInv, "try-error")){
      avInfInv <- try(solve(avInf + diag(tolParInv,ncol(avInf))), silent = TRUE)
    }
    ##########################
    # 5) get 1st derivatives from MME-version (correct)
    # PAPER FORMULA (Lee and Van der Werf, 2006)
    # dL/ds2u = -0.5 [(Nu/s2u) - (tr(AiCuu)/s4u) -  (e/s2e)'(Zu/s2u)]
    # dL/ds2e = -0.5 [((Nr-Nb)/s2e) - [(Nu - (tr(AiCuu)/s2u))*(1/s2e)] - ... - (e/s2e)'(e/s2e)]
    # 
    # PAPER FORMULA (Jensen and Madsen, 1997)
    # dL/ds2u = (q.i * lambda) - (lambda * (T + S) * lambda)  Eq. 18
    # dL/ds2e = tr(Rij*Ri) - tr(Ci*W'*Ri*Rij*Ri*W) - (e'*Ri*Rij*Ri*e)
    ###########################
    
    # invert the coefficient matrix
    Ci=solve(MChol[-ncol(MChol),-ncol(MChol)]) # this is the cholesky decomposition of C
    Ci = (Ci)%*%t(Ci) # multiply by its transpose since the decomposition is only one of the triangular parts
    Ci <- as.matrix(Ci)
    
    emInfInvList <- emInfList <- list() # store EM information matrices and its inverses
    # calculate first derivatives and take advantage to get emInf and emInfInv
    dLu = list() # store first derivatives
    
    if(is.null(Z)){ # if random effects exist
      iR=0
    }else{
      for(iR in 1:length(zsAva)){ # iR=1 # for each random effect
        SigmaInv = lambda[[iR]] # get the inverse of the variance-cov components for this random effect
        omega<-omega2<-traces <- matrix(0,nrow = nrow(SigmaInv),ncol=ncol(SigmaInv))
        # calculate T or tr(AiCii) Eq. 18.1 of Jensen
        for(k in 1:nrow(SigmaInv)){
          for(l in 1:nrow(SigmaInv)){
            if(thetaC[[iR]][k,l]>0){
              usedPartitionK = partitions[[iR]][k,1]:partitions[[iR]][k,2]
              usedPartitionL = partitions[[iR]][l,1]:partitions[[iR]][l,2]
              trAiCuu = sum(diag((Ai[[iR]])%*%Ci[usedPartitionK,usedPartitionL]))
              traces[k,l] = trAiCuu
            }else{
              traces[k,l] = 0
            }
          }
        }
        traces[lower.tri(traces)] <- t(traces)[lower.tri(traces)]
        ## first derivatives = dL/ds2u = (q.i * lambda) - (lambda * (T + S) * lambda)    where S=UAiU and we use U.lambda
        dLuProv = ((Nus[iR])*SigmaInv) -  (t(U.SinvList[[iR]])%*%Ai[[iR]]%*%U.SinvList[[iR]]) - (SigmaInv%*%traces%*%SigmaInv)
        ## althernative EM update
        # current(theta)   -   update(delta)  but we need to decompose the update(delta) = Iem * vech(dLu/ds2u) , Iem is then of dimensions equal to vech(dLu/ds2u)
        # theta[[iR]] - (theta[[iR]]%*%dLuProv%*%theta[[iR]])/Nus[iR]    Eq.34 
        thetaUnlisted <- theta[[iR]][which(thetaC[[iR]] > 0)]
        thetaUnlistedMat = diag(thetaUnlisted,length(thetaUnlisted),length(thetaUnlisted))
        # emInfInvProvExt <- ( thetaUnlisted %*% t(thetaUnlisted) )/Nus[iR]
        emInfInvProvExt <- thetaUnlistedMat %*% t(thetaUnlistedMat)/Nus[iR]
        emInfInvList[[iR]] <- emInfInvProvExt
        emInfList[[iR]] <- solve(emInfInvList[[iR]] ) 
        dLu[[iR]] = dLuProv[which(thetaC[[iR]]>0)] # vech(dLu/ds2u) or half-vectoization
      }
    }
    
    dLe <- numeric() # first derivatives for error variance components Eq. 19 of Jensen
    for(iS in 1:length(S)){ # Rij <- S[[iS]]%*%Ri
      dLe[iS] = (
        (sum(diag(S[[iS]]%*%Ri)) - sum(diag(Ci%*%t(W)%*%Ri%*%S[[iS]]%*%Ri%*%W)) ) - 
          ( t(e)%*%Ri%*%S[[iS]]%*%Ri%*%e )
      )
    }
    dLu[[iR+1]] <- dLe
    dLeM <- vectorToList(thetaC[[ep]],dLe)
    
    thetaRUnlisted <- theta[[iR+1]][which(thetaC[[iR+1]] > 0)]
    thetaRUnlistedMat = diag(thetaRUnlisted,length(thetaRUnlisted),length(thetaRUnlisted))
    emInfInvProvResExt <- thetaRUnlistedMat %*% t(thetaRUnlistedMat)/Nr
    emInfInvList[[iR+1]] <-emInfInvProvResExt
    emInfList[[iR+1]] <- solve(emInfInvList[[iR+1]])
    
    dLu = unlist(dLu); # first derivatives as vector
    emInf = do.call(adiag1,emInfList) # big EM information matrix
    # cons = do.call(adiag1,thetaC) # big constraint matrix
    ###########################
    # 6) update the variance paramters (CORRECT)
    # PAPER FORMULA (Lee and Van der Werf, 2006)
    # theta.n+1 = theta.n + (AInfi * dL/ds2) 
    ###########################
    
    # unlist parameters and constraints
    thetaUnlisted <- unlistThetaWithThetaC(theta,thetaC)
    thetaCUnlisted <- unlistThetaWithThetaC(thetaC,thetaC)
    
    weightInfMatCurrentIter=weightInfMat[iIter] #[iIter+stepsEM]
    # Joint information matrix and update
    #                   AVERAGE INFORMATION                       +     EXPECTATION MAXIMIZATION
    InfMat = ((diag(ncol(avInf))-weightInfEMl[[iIter]]) %*% avInf ) + (weightInfEMl[[iIter]] %*% emInf) # combined information matrix EM + avInf
    InfMatInv <- solve(InfMat) # inverse of the information matrix
    change = (InfMatInv*weightInfMatCurrentIter)%*%as.matrix(dLu) # delta = I- * dLu/dLx
    # print(length(which(thetaCUnlisted == 3)))
    expectedNewTheta = thetaUnlisted - change
    
    ####################
    ## 1) apply constraints
    
    for(ivc in 1:length(thetaCUnlisted)){
      if(thetaCUnlisted[ivc] == 1){
        if(expectedNewTheta[ivc] < 1e-10){
          cat("constraining to small value")
          expectedNewTheta[ivc]=1e-10
          toBoundary[iIter,ivc]=1 # keep track of those VC that are going to the boundaries
          ## change to fixed
          sumToBoundary <- apply(toBoundary,2,sum)
          toForce <- which(sumToBoundary >= 3)
          if(length(toForce) > 0){
            thetaCUnlisted[toForce]=3
          }
        }}
      if(thetaCUnlisted[ivc] == 3){
        thetaUnlistedPlusAddScaleParam <- c(expectedNewTheta,addScaleParam)
        # theta.i            = scaleParameter.selected      *  Theta
        expectedNewTheta[ivc]=thetaF[ivc,] %*%  thetaUnlistedPlusAddScaleParam
      }
    }
    
    ####################
    ## 2) if there's fixed vc use a different update
    if(length(which(thetaCUnlisted == 3)) > 0){
      cat("updates using constraints")
      InfMat.uu = InfMat[which(thetaCUnlisted != 3),which(thetaCUnlisted != 3)]
      InfMat.ff = InfMat[which(thetaCUnlisted == 3),which(thetaCUnlisted == 3)]
      dLu.uu =  dLu[which(thetaCUnlisted != 3)]
      dLu.ff =  dLu[which(thetaCUnlisted == 3)]
      InfMat.uf = InfMat[which(thetaCUnlisted != 3),which(thetaCUnlisted == 3)]
      InfMatInv.uu <- ginv(InfMat.uu) # inverse of the information matrix
      InfMatInv.ff <- ginv(InfMat.ff) # inverse of the information matrix
      change.ff = InfMatInv.ff %*% dLu.ff
      change.uu = InfMatInv.uu %*% dLu.uu #(dLu.uu - (InfMat.uf%*%change.ff) )
      change[which(thetaCUnlisted != 3)]= change.uu
      change[which(thetaCUnlisted == 3)]= 0#change.ff
      ### new values for variance components theta.i+1 = theta.i + delta
      expectedNewTheta = thetaUnlisted - change
    }
    
    ####################
    # 3) quantify change and norm of delta ||D||
    if(iIter > 1){
      percChangeMat[,iIter] <- (change/change0)*100
      changeCheck <- which(abs(percChangeMat[,iIter]) > 15000)
      changeCheck <- setdiff(changeCheck,which(toBoundary[iIter,] > 0))
      # if(length(changeCheck) > 0){ # if there's vc that changed drastically decrease the update for that VC
      #   cat(red("Reducing Updates by 0.1"))
      #   change[changeCheck] <- change[changeCheck] * .1
      #   expectedNewTheta = thetaUnlisted - change
      # }
      change0 <- change
    }else{
      change0 <- change
      percChangeMat[,iIter] <- 0
    }
    # cat(green(paste("change",paste(round(percChangeMat[,iIter],0), collapse = "%, "))),"\n")
    #########################################################
    ### 4) if not positive definite change to full EM update or add a value to the diagonal
    pdCheck <- vectorToList(thetaC, expectedNewTheta)
    pdCheck <- do.call(adiag1,pdCheck)
    pdCheck <- eigen(pdCheck)$values
    pdCheck <- which(pdCheck < 0)
    if(length(pdCheck) > 0){ # if new matrix is not positive definite use the EM update
      # cat(red("Updated VC matrix is not positive definite, adding a small value to diagonal"))
      # InfMatInv <- solve(InfMat + diag(1e-5, ncol(InfMat), ncol(InfMat) )) # add a value to the diagonal
      cat(red("Updated VC matrix is not positive definite, changing to an EM step"))
      InfMat = (0.5 * avInf ) + (0.5 * emInf) # combined information matrix EM + avInf # change to EM update
      
      if(length(which(thetaCUnlisted == 3)) > 0){
        cat("updates using constraints in an EM step")
        InfMat.uu = InfMat[which(thetaCUnlisted != 3),which(thetaCUnlisted != 3)]
        InfMat.ff = InfMat[which(thetaCUnlisted == 3),which(thetaCUnlisted == 3)]
        dLu.uu =  dLu[which(thetaCUnlisted != 3)]
        dLu.ff =  dLu[which(thetaCUnlisted == 3)]
        InfMat.uf = InfMat[which(thetaCUnlisted != 3),which(thetaCUnlisted == 3)]
        InfMatInv.uu <- ginv(InfMat.uu) # inverse of the information matrix
        InfMatInv.ff <- ginv(InfMat.ff) # inverse of the information matrix
        change.ff = InfMatInv.ff %*% dLu.ff
        change.uu = InfMatInv.uu %*% dLu.uu #(dLu.uu - (InfMat.uf%*%change.ff) )
        change[which(thetaCUnlisted != 3)]= change.uu
        change[which(thetaCUnlisted == 3)]= 0#change.ff
        ### new values for variance components theta.i+1 = theta.i + delta
        expectedNewTheta = thetaUnlisted - change
      }else{
        InfMatInv = ginv(InfMat)
        expectedNewTheta = thetaUnlisted - ( InfMatInv%*% as.matrix(dLu) )
      }
    }
    
    ## stoping criterias
    changeNorm[1,iIter]=norm(as.matrix(change[which(thetaCUnlisted != 3)])) # stopping criteria 1
    changeNorm[2,iIter]=norm(as.matrix(dLu[which(thetaCUnlisted != 3)])) # stopping criteria 2
    changeNorm[3,iIter]= norm(as.matrix( (diag(InfMatInv[which(thetaCUnlisted != 3),which(thetaCUnlisted != 3)]))/sqrt(length(dLu[which(thetaCUnlisted != 3)])) * dLu[which(thetaCUnlisted != 3)] )) # stopping criteria 3
    
    ####################
    # bring VC from vector to matrices
    expectedNewThetaList = vectorToList(thetaC, expectedNewTheta)
    for(i in 1:length(expectedNewThetaList)){
      expectedNewThetaList[[i]] <- as.matrix(nearPD(expectedNewThetaList[[i]])$mat)
    }
    # print(expectedNewThetaList)
    # update covariance parameters and save for monitoring
    theta <- expectedNewThetaList
    finalThetaList[[iIter]] <- expectedNewTheta
    
    ####################
    # stopping criteria
    if(iIter > 1){
      delta_llik = llik[iIter] - llik[iIter-1]
      ## stopping criteria 0
      # if( (  (delta_llik < tolParConvNorm)) | (iIter == nIters)  ){ 
      #   break 
      # }
      # # stopping criteria 1
      # if( (  (changeNorm[1,iIter] < 0.05)) | (iIter == nIters)  ){
      #   break
      # }
      # stopping criteria 3
      if( (  ( (changeNorm[3,iIter] < tolParConvNorm)) & (changeNorm[1,iIter] < tolParConvNorm) ) | (iIter == nIters) | (delta_llik < tolParConvLL) ){
        break
      }
      
    }
  } # end of AI
  monitor <- do.call(cbind,finalThetaList)
  # plot(llik)
  rownames(monitor) <- c(paste0("vc",1:(nrow(monitor))))
  return(list(sigma=theta, sigmavector=expectedNewTheta, yPy=yPy, llik=llik, b=b, u=u, Wu=Wu,  
              avInf=avInf, emInf=emInf, Ci=Ci, dLu=dLu, 
              M=M, C=C, monitor=monitor, normChange3=normChange3,
              logDetA=logDetA, Nus=Nus, logDetC=logDetC, yPy=yPy, Wu=Wu, change=change,
              expectedNewTheta=expectedNewTheta, expectedNewThetaList=expectedNewThetaList,
              percChange=percChangeMat, changeNorm=changeNorm,thetaCUnlisted=thetaCUnlisted,
              sumToBoundary=sumToBoundary))
}

unlistThetaWithThetaC <- function(x, xc){
  thetaUnlisted= unlist(lapply(as.list(1:length(x)), function(v){
    return(x[[v]][which(xc[[v]]>0)])
  }))
  return(thetaUnlisted)
}


pmonitor <- function(object,...){
  x <- object$monitor
  mmin <- min(x)
  mmax <- max(x)
  plot(x[1,], type="l", ylim=c(mmin,mmax))
  if(nrow(x) > 1){
    for(i in 2:nrow(x)){
      par(new=TRUE)
      plot(x[i,], col=i, ylim=c(mmin,mmax), ylab="",xlab="", type="l",...)
    }
  }
  legend("topright",legend = rownames(x), col=1:(nrow(x)), bty="n", lty=1)
}

vectorToList <- function(constraint, parameters){
  parametersList <- constraint
  start <- 1
  for(o in 1:length(constraint)){
    last <- start+length(which(constraint[[o]]>0))-1
    provParam <- parametersList[[o]] # provisional covariance matrix
    provParam[which(constraint[[o]]>0)] <- parameters[start:last] # fill covariance matrix with new parameters
    provParam[lower.tri(provParam)] <- t(provParam)[lower.tri(provParam)] # fill lower triangular
    parametersList[[o]] <- provParam # return new covariance matrix to original object
    start <- last+1 # update start
  }
  return(parametersList)
}

mmeFormation <- function(X,Z=NULL,Ri=NULL,y){
  # X is a matrix
  # Z is a list of lists
  # y is a vector
  # Zlist <- unlist(Z) # put Z in a list of one level
  if(is.null(Z)){
    XZlist <- c(list(X))
    XZylist <- c(list(X),list(matrix(y))) 
  }else{
    # Zlist <- lapply(rapply(Z, enquote, how="unlist"), eval)
    XZlist <- c(list(X),Z)
    XZylist <- c(list(X),Z,list(matrix(y))) 
  }
  
  W <- as.matrix(do.call(cbind, XZlist))
  
  # Zlist <- NULL # delete to avoid memory consumption
  ncolM <- unlist(lapply(XZylist, ncol))
  end <- numeric()
  for(u in 1:length(ncolM)){end[u] <- sum(ncolM[1:u])}
  start <- end-ncolM+1
  ncolMsum <- sum(ncolM)
  M <- matrix(NA,ncolMsum,ncolMsum)
  #--->
  # -->
  #  ->
  for(i in 1:length(XZylist)){ # for each row of the MME
    useRow <- start[i]:end[i]
    for(j in i:length(XZylist)){ # for each col of the MME
      useCol <- start[j]:end[j]
      if(is.null(Ri)){
        M[useRow,useCol] <- as.matrix(t(XZylist[[i]])%*%XZylist[[j]])
        M[useCol,useRow] <- t(M[useRow,useCol])
      }else{
        M[useRow,useCol] <- as.matrix(t(XZylist[[i]])%*%Ri%*%XZylist[[j]])
        M[useCol,useRow] <- t(M[useRow,useCol])
      }
    }
  }
  return(list(M=M, W=W, start=start, end=end))
}

image2 <- function(x){
  print(image(as(x, Class = "sparseMatrix")))
}

# emInfInvProv <- (((theta[[iR]]%*%dLuProv%*%theta[[iR]]))%*%solve(dLuProv))/Nus[iR] # real em information matrix
# print(dim(emInfInvProvExt))
# emInfInvProv <- (((theta[[iR]]%*%t(theta[[iR]]))))/Nus[iR] # real em information matrix
# 
# indicator <- as.data.frame(which(thetaC[[iR]] > 0, arr.ind = TRUE));  # indicator df to tranform emInfInvProv into same dimensions of the AI matrix
# emInfInvProvExt <- matrix(0, Nvc2[iR], Nvc2[iR]) # extended EM inf matrix
# for(ii in 1:nrow(indicator)){ # for each variance-covariance component
#   SeconDer <- indicator[ii,1] # which rows of the emInfInvProv (information matrix) to take and put into the ii'th row of emInfInvProvExt
#   firstDer <- indicator[ii,2] # which columns of first derivatives will be multiplied by
#   pickCol <- sort(unique(c(which(indicator[,1]==firstDer),which(indicator[,2]==firstDer)))) # vc-cov components that should be multiplied
#   pickColInemInfInvProv <- which(thetaCs[[iR]][SeconDer,] > 0) # only take elements in the information matrix that match elements in theta to be estimated
#   emInfInvProvExt[ii,pickCol] = emInfInvProv[SeconDer,pickColInemInfInvProv] # add elements ofemInfInvProv into emInfInvProvExt
# }

# emInfInvProvResExt <- (((theta[[iR+1]]%*%dLeM%*%theta[[iR+1]]))%*%solve(dLeM))/Nr # real em.r information matrix
# 
# indicatorR <- as.data.frame(which(thetaC[[ep]] > 0, arr.ind = TRUE)); # indicator df to tranform emInfInvProv into same dimensions of the AI matrix
# emInfInvProvResExt <- matrix(0, Nvc2[ep], Nvc2[ep]) # extended EM.r inf matrix
# for(ii in 1:nrow(indicatorR)){ # # for each residual variance-covariance component
#   SeconDer <- indicatorR[ii,1] # which rows of the emInfInvProv (information matrix) to take and put into the ii'th row of emInfInvProvExt
#   firstDer <- indicatorR[ii,2] # which columns of first derivatives will be multiplied by
#   pickCol <- sort(unique(c(which(indicatorR[,1]==firstDer),which(indicatorR[,2]==firstDer))))  # vc-cov components that should be multiplied
#   pickColInemInfInvProvRes <- which(thetaCs[[ep]][SeconDer,] > 0) # only take elements in the information matrix that match elements in theta to be estimated
#   emInfInvProvResExt[ii,pickCol] = emInfInvProvRes[SeconDer,pickColInemInfInvProvRes] # add elements ofemInfInvProv into emInfInvProvExt
# }
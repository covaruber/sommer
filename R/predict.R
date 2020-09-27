#### =========== ######
## PREDICT FUNCTION #
#### =========== ######

"predict.mmer" <- function(object,classify=NULL,...){
  oto <- object$terms
  for(u in 1:length(object$terms)){
    prov <- object$terms[[u]]
    if(length(prov) > 0){
      for(v in 1:length(prov)){
        object$terms[[u]][[v]] <- unlist(strsplit(prov[[v]],":"))
      }
    }
  }
  ##################################################
  # step 0. find all variables used in the modeling
  allTermsUsed <- unique(c(unlist(object$terms$fixed), unlist(object$terms$random)))
  allTermsUsed<- allTermsUsed[which(allTermsUsed!= "1")]
  allTermsUsed <- unique(unlist(strsplit(allTermsUsed,":")))
  
  toAgg <- unique(unlist(strsplit(classify,":")))
  ignored <- setdiff(allTermsUsed,toAgg)
  
  levelsOfTerms <- apply(data.frame(toAgg),1,function(x){unique(object$dataOriginal[,x])})
  DTX <- expand.grid(levelsOfTerms); colnames(DTX) <- toAgg
  DTX <-merge(object$dataOriginal[,c(colnames(DTX),ignored,object$terms$response[[1]])], DTX, all.y = TRUE)
  
  # DTX <-merge(DTX, object$dataOriginal[,c(colnames(DTX),ignored)], all.x = TRUE)
  
  if(length(object$terms$response[[1]]) < 2){
    YY = data.frame(DTX[,object$terms$response[[1]]]); colnames(YY) <- object$terms$response[[1]]
  }else{YY = DTX[,object$terms$response[[1]]]}
  
  DTX[,object$terms$response[[1]]] <- apply(YY,2,imputev)
  if(length(ignored) == 1){
    # v=which(colnames(DTX) == ignored)
    DTX[,ignored] <- imputev(DTX[,ignored])
  }else{
    DTX[,ignored] <- apply(DTX[,ignored],2,imputev)
  }
  
  ##################################################
  # step 1 create all models
  # 1. extended data model (just get matrices)
  # 2. original model without reshaped output
  # 3. original model (just get matrices)
  if(is.null(object$call$random)){
    modelForMatrices <- mmer(fixed=object$call$fixed,
                             # random=object$call$random,
                             rcov=object$call$rcov, 
                             data=DTX, return.param = TRUE,#reshape.results=TRUE,
                             na.method.Y = object$call$na.method.Y,
                             na.method.X = object$call$na.method.X,...)
    originalModel <- mmer(fixed=object$call$fixed,
                          # random=object$call$random,
                          rcov=object$call$rcov,
                          data=object$dataOriginal, return.param = FALSE,reshape.output =FALSE,
                          init = object$sigma_scaled, constraints = object$constraints,
                          na.method.Y = object$call$na.method.Y, 
                          na.method.X = object$call$na.method.X,...)
    originalModelForMatricesSE <- mmer(fixed=object$call$fixed,
                                       # random=object$call$random,
                                       rcov=object$call$rcov,
                                       data=object$dataOriginal, return.param = TRUE,#reshape.output =FALSE,
                                       init = object$sigma_scaled, constraints = object$constraints,
                                       na.method.Y = object$call$na.method.Y,
                                       na.method.X = object$call$na.method.X,...)
  }else{
    modelForMatrices <- mmer(fixed=object$call$fixed,
                             random=object$call$random,
                             rcov=object$call$rcov, 
                             data=DTX, return.param = TRUE,#reshape.results=TRUE,
                             na.method.Y = object$call$na.method.Y, 
                             na.method.X = object$call$na.method.X, ...)
    originalModel <- mmer(fixed=object$call$fixed,
                          random=object$call$random,
                          rcov=object$call$rcov,
                          data=object$dataOriginal, return.param = FALSE,reshape.output =FALSE,
                          init = object$sigma_scaled, constraints = object$constraints,
                          na.method.Y = object$call$na.method.Y,
                          na.method.X = object$call$na.method.X,...)
    originalModelForMatricesSE <- mmer(fixed=object$call$fixed,
                          random=object$call$random,
                          rcov=object$call$rcov,
                          data=object$dataOriginal, return.param = TRUE,#reshape.output =FALSE,
                          init = object$sigma_scaled, constraints = object$constraints,
                          na.method.Y = object$call$na.method.Y,
                          na.method.X = object$call$na.method.X,...)
  }
  modelForMatrices$U <- originalModel$U
  modelForMatrices$PevU <- originalModel$PevU
  modelForMatrices$VarU <- originalModel$VarU
  modelForMatrices$Beta <- originalModel$Beta
  modelForMatrices$VarBeta <- originalModel$VarBeta
  modelForMatrices$Vi <- originalModel$Vi
  modelForMatrices$P <- originalModel$P
  modelForMatrices$sigma <- originalModel$sigma
  ##################################################
  # calculate Xb and get X multivariate matrices from extended and original model
  ys <- object$terms$response[[1]]
  nt <- length(ys) # number of traits
  TT <- diag(nt) # diagonal matrix
  # which fixed effects to include
  fToUse <- list()
  for(i in 1:length(toAgg)){
    fToUse[[i]]<-grep(toAgg[i],object$terms$fixed[[1]])
  }
  fToUse = sort(c(1,unique(unlist(fToUse))))
  X <- do.call(cbind,modelForMatrices$X[fToUse]) # build X cbinding the ones required
  # find the betas to use
  ncolsX <- unlist(lapply(modelForMatrices$X,ncol))
  start=1; betas0 <- list()
  for(i in 1:length(ncolsX)){
    betas0[[i]] <- start:(start+(ncolsX[i]*nt)-1)
    start= max(betas0[[i]])+1
  }
  # new X extended
  Xm.extended <- kronecker(X,TT)
  Xb=Xm.extended%*%modelForMatrices$Beta[unlist(betas0[fToUse]),1] # calculate Xb
  # X'ViX
  XtViX = modelForMatrices$VarBeta[unlist(betas0[fToUse]),unlist(betas0[fToUse])]
  # build X multivariate from original model
  Xo <- do.call(cbind,originalModelForMatricesSE$X[fToUse])
  Xm.original <- kronecker(Xo,TT)
  ##################################################
  # calculate Zu
  if(!is.null(object$call$random)){
    nre <- length(object$terms$random) # number of random effects
    reUsed0 <- list()
    for(i in 1:length(toAgg)){
      reUsed <- numeric()
      for(j in 1:nre){
        result <- grep(toAgg[i],object$terms$random[[j]])
        if(length(result) >0){reUsed[j] =1}else{reUsed[j] =0}
      }
      reUsed0[[i]] <- reUsed
    }; 
    reUsed0 <- Reduce("+",reUsed0); reUsed0[which(reUsed0 > 0)]=reUsed0[which(reUsed0 > 0)]/reUsed0[which(reUsed0 > 0)]
    # identify which effects estimated correspond to such random terms
    zToUse <- list(); start=1
    for(i in 1:nre){
      zToUse[[i]] <-  start:(start+object$termsN$random[i]-1)
      start <- max(zToUse[[i]])+1
    }
    zToUse <- unique(unlist(zToUse[which(reUsed0 == 1)]))
    
    nz <- length(modelForMatrices$Z)
    Zu <- vector(mode = "list", length = nre) # list for Zu
    for(ir in 1:nz){ # for each random effect
      Z <- modelForMatrices$Z[[ir]] # provisional Z
      Zu[[ir]] <- kronecker(Z,TT) %*% modelForMatrices$U[[ir]] # calculate Zu
    }
    ###########################################
    ## SE step
    # cov(b,u - u.hat), xvxi, pev (C12,C11,C22)
    ## build the multivariate K and Z from the original odel to build 
    # ZKfv, Xo, cov.b.pev, pev
    if(!is.null(zToUse)){ # if not only there's random terms but they are relevant for predict
      VarKl <- lapply(as.list((zToUse)),function(x){
        kronecker(as.matrix(originalModelForMatricesSE$K[[x]]),originalModel$sigma[,,x])
      }); VarK <- do.call(adiag1,VarKl)
      
      Zl <- lapply(as.list((zToUse)),function(x){
        kronecker(as.matrix(t(originalModelForMatricesSE$Z[[x]])),TT)
      }); tZm <- do.call(rbind,Zl)
      
      ZKfv <- VarK %*% tZm # as many rows as obs, as many cols as levels
      
      cov.b.pev <- 0 - t(originalModel$P %*% Xm.original) %*% originalModel$Vi %*% t(ZKfv)
      pev <- do.call(adiag1,originalModel$PevU[zToUse])
      # (185 x 3)' (185 x 185) = (3 x 185)  (185 x 185) (levs164 x obs185)' = (3 x 164)
      # bring the design matrices for extended model to get standard errors for the extended model
      Znew <- lapply(as.list((zToUse)),function(x){
        kronecker(as.matrix(t(modelForMatrices$Z[[x]])),TT)
      }); Znew <- do.call(rbind,Znew)
      tZnew <- t(Znew)
    }
  }
  ##################################################
  # add them up
  if(!is.null(object$call$random) & !is.null(zToUse)){
    y.hat <- Xb + Reduce("+",Zu[zToUse]) # y.hat = Xb + Zu.1 + ... + Zu.n
    standard.errors <- sqrt(rowSums((Xm.extended%*%XtViX)*Xm.extended) +
                              rowSums(2*(Xm.extended%*%cov.b.pev)*tZnew) +
                              rowSums((tZnew%*%pev)*tZnew))
  }else{
    y.hat <- Xb # y.hat = Xb 
    standard.errors <- sqrt(rowSums((Xm.extended%*%XtViX)*Xm.extended))
  }
  ##################################################
  # add y.hat to the grid
  
  if(length(allTermsUsed) == 1){
    DTX2 <- data.frame(DTX[,-which(colnames(DTX) %in% ys)])
    colnames(DTX2) <- allTermsUsed
  }else{
    DTX2 <- DTX[,-which(colnames(DTX) %in% ys)]
  }
  DTX2L <- rep(list(DTX2),nt)
  for(it in 1:nt){DTX2L[[it]]$trait <- ys[it]}
  DTX2 <- do.call(rbind, DTX2L)
  realOrder <- list()
  for(i in 1:nt){
    realOrder[[i]] <- seq(i,length(y.hat[,1]),nt)
  }
  realOrder <- unlist(realOrder)
  DTX2$predicted.value <- y.hat[realOrder,1]
  DTX2$standard.error <- standard.errors[realOrder]
  DTX2$status <- "Estimable"
  DTX2$status[which(is.nan(DTX2$standard.error))] <- "Issues"
  DTX2$standard.error[which(is.nan(DTX2$standard.error))] <- 0
  
  ##################################################
  # aggregate to desired shape
  
  myForm <- paste0("predicted.value~",paste(toAgg, collapse = "+"), "+trait")
  pvals <- aggregate(as.formula(myForm), FUN=mean, data=DTX2)
  myFormSE <- paste0("standard.error~",paste(toAgg, collapse = "+"), "+trait")
  pvalsSE <- aggregate(as.formula(myFormSE), FUN=mean, data=DTX2)
  # pvalsSE <- aggregate(as.formula(myForm), FUN=function(x){1.96 * (sd(x)/sqrt(length(x)))}, data=DTX2)
  colnames(pvalsSE)[ncol(pvalsSE)] <- "standard.error"
  pvals <- merge(pvals,pvalsSE, by=c("trait",toAgg))
  
  toreturn2 <- list(pvals=pvals, 
                    FE.used = unique(unlist(object$terms$fixed[[1]][fToUse])), 
                    RE.used=unique(unlist(oto$random[which(reUsed0==1)])),
                    Ignored= setdiff(allTermsUsed,toAgg),
                    model=modelForMatrices
                    #EE.used=unique(unlist(object$terms$rcov))
  )
  attr(toreturn2, "class")<-c("predict.mmer", "list")
  
  return(toreturn2)
}


"print.predict.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  aver <- paste(x$Ignored,collapse = ", ")
  fused <- paste(x$FE.used,collapse = ", ")
  rused <- paste(x$RE.used,collapse = ", ")
  cat(blue(paste("\n  The predictions are obtained by averaging across the hypertable
                 calculated from model terms constructed solely from factors in
                 the averaging and classify sets.
                 - The simple averaging set:",aver,"\n",
                 "- The fixed effects included:",fused,"\n",
                 "- The random effects included:",rused,"\n\n")
  ))
  head(x$pvals,...)
  # head(print.predict.mmer(pp))
  # print((x$predictions))
}

"head.predict.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  aver <- paste(x$Ignored,collapse = ", ")
  fused <- paste(x$FE.used,collapse = ", ")
  rused <- paste(x$RE.used,collapse = ", ")
  cat(blue(paste("\n  The predictions are obtained by averaging across the hypertable
                 calculated from model terms constructed solely from factors in
                 the averaging and classify sets.
                 - The simple averaging set:",aver,"\n",
                 "- The fixed effects included:",fused,"\n",
                 "- The random effects included:",rused,"\n\n")
  ))
  head(x$pvals, ...)
  # head(print.predict.mmer(pp))
  # print((x$predictions))
}


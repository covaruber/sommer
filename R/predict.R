#### =========== ######
## PREDICT FUNCTION #
#### =========== ######
# include is used for aggregating
# averaged is used to be included in the prediction
# ignored is not used included in the prediction
"predict.mmer" <- function(object,classify=NULL,hypertable=NULL,...){
  
  classify<- unique(unlist(strsplit(classify,":")))

  if(is.null(classify)){
    stop("Please provide the classify argument. For fitted values use the fitted() function.",call. = FALSE)
  }

  oto <- oto2 <- object$terms
  oto2$fixed[[1]] <- setdiff(oto2$fixed[[1]],c("1","-1"))
  oto2$fixed <- lapply(oto2$fixed,function(x){paste(x,collapse = ":")})
  oto2$random <- lapply(oto2$random,function(x){paste(x,collapse = ":")})
  
  for(u in 3:length(object$terms)){ # change random terms to split by ":"
    prov <- object$terms[[u]]
    if(length(prov) > 0){
      for(v in 1:length(prov)){
        object$terms[[u]][[v]] <- unlist(strsplit(prov[[v]],":"))
      }
    }
  }
  
  include <- setdiff(unique(c(unlist(object$terms$fixed),unlist(object$terms$random))),c("1","-1"))
  
  if(!is.null(hypertable)){ # if user provides a hypertable, use the customization instead
    include <- setdiff(hypertable[which(hypertable$include),"termHL"],c("1","-1"))
    # print(include)
  }
  
  ##################################################
  # step 0. find all variables used in the modeling
  allTermsUsed <- unique(c(unlist(object$terms$fixed), unlist(object$terms$random)))
  allTermsUsed<- allTermsUsed[which(allTermsUsed!= "1")]
  allTermsUsed<- allTermsUsed[which(allTermsUsed!= "-1")]
  allTermsUsed <- unique(unlist(strsplit(allTermsUsed,":")))
  
  toAgg <- unique(unlist(strsplit(include,":")))
  ignored <- setdiff(allTermsUsed,toAgg)
  # print(ignored)
  
  levelsOfTerms <- apply(data.frame(toAgg),1,function(x){unique(object$dataOriginal[,x])})
  DTX <- expand.grid(levelsOfTerms); 
  colnames(DTX) <- toAgg
  # toMerge <- object$dataOriginal[,c(colnames(DTX),ignored,object$terms$response[[1]])]
  
  toMerge <- unique(object$dataOriginal[,c(colnames(DTX),ignored)])
  if(!is.data.frame(toMerge)){
    toMerge <- data.frame(toMerge); colnames(toMerge) <- colnames(DTX)
  }
  toMerge[,object$terms$response[[1]]] <- 1
  # print(toMerge)
  # print(DTX)
  DTX <-merge(toMerge, DTX, all.y = TRUE) # by=intersect(colnames(DTX),colnames(toMerge)),
  # print(str(DTX))
  # DTX <-merge(DTX, object$dataOriginal[,c(colnames(DTX),ignored)], all.x = TRUE)
  
  if(length(object$terms$response[[1]]) < 2){
    YY = data.frame(DTX[,object$terms$response[[1]]]); colnames(YY) <- object$terms$response[[1]]
  }else{YY = DTX[,object$terms$response[[1]]]}
  
  DTX[,object$terms$response[[1]]] <- apply(YY,2,imputev)
  if(length(ignored) > 0){ # if there's ignored columns
    if(length(ignored) == 1){
      # v=which(colnames(DTX) == ignored)
      DTX[,ignored] <- imputev(DTX[,ignored])
    }else{
      for(o in 1:length(ignored)){
        DTX[,ignored[o]] <- imputev(DTX[,ignored[o]])
      }
    }
  }
  # print(DTX)
  # print(dim(DTX))
  # print(str(DTX))
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
  # print(object$terms$fixed[[1]])
  fToUse <- list()
  for(i in 1:length(toAgg)){
    fToUse[[i]]<-grep(toAgg[i],object$terms$fixed[[1]])
  }
  fToUse = sort(c(1,unique(unlist(fToUse))))
  # print(fToUse)
  # print(length(modelForMatrices$X))
  if(!is.null(hypertable)){ # if user provides a hypertable, use the customization instead
    fToUse <- which(hypertable$type == "fixed" & hypertable$include == TRUE)
    fNotToUse <- which(hypertable$type == "fixed" & hypertable$ignored == TRUE)
    if(length(fNotToUse) > 0){fToUse <- setdiff(fToUse,fNotToUse)}
  }
  
  if(length(fToUse) > 0){ # if there's fixed effects to add
    X <- do.call(cbind,modelForMatrices$X[fToUse]) # build X cbinding the ones required
    # print(head(X))
    # find the betas to use
    ncolsX <- unlist(lapply(modelForMatrices$X,ncol))
    # print(ncolsX)
    start=1; betas0 <- list()
    for(i in 1:length(ncolsX)){
      if(ncolsX[i] > 0){
        betas0[[i]] <- start:(start+(ncolsX[i]*nt)-1)
        start= max(betas0[[i]])+1
      }
    }
    # print(betas0)
    # new X extended
    X.mv.extended <- kronecker(X,TT)
    Xb=X.mv.extended%*%modelForMatrices$Beta[unlist(betas0[fToUse]),1] # calculate Xb
    # X'ViX
    XtViX = modelForMatrices$VarBeta[unlist(betas0[fToUse]),unlist(betas0[fToUse])]
    # build X multivariate from original model
    Xo <- do.call(cbind,originalModelForMatricesSE$X[fToUse])
    X.mv.original <- kronecker(Xo,TT)
  }
  ##################################################
  # calculate Zu
  pev=NULL
  cov.b.pev=NULL
  zToUse=NULL
  if(!is.null(object$call$random)){
    nre <- length(object$terms$random) # number of random effects
    # identify which random terms in the model should be added
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
    # print(reUsed0)
    
    # identify which effects estimated correspond to each random term in the model
    zToUse <- list(); start=1
    for(i in 1:nre){
      zToUse[[i]] <-  start:(start+object$termsN$random[i]-1)
      start <- max(zToUse[[i]])+1
    }
    zToUse <- unique(unlist(zToUse[which(reUsed0 == 1)]))
   
    
    if(!is.null(hypertable)){ # if user provides a hypertable, use the customization instead
      nFixed <- max(which(hypertable$type == "fixed"))
      # which ones the user wants to force to be included
      rForce <- which(hypertable$type == "random" & hypertable$include == TRUE)
      if(length(rForce) > 0){
        rForce <- rForce - nFixed
        zToUse <- rForce
      }
      # which ones the user wants to force to be excluded
      rNotForce <- which(hypertable$type == "random" & hypertable$ignored == TRUE)
      if(length(rNotForce) > 0){
        rNotForce <- rNotForce - nFixed
        zToUse <- setdiff(zToUse,rNotForce)
      }
    }
    
    if(length(zToUse) == 0){zToUse=NULL}
    # print(zToUse)
    nz <- length(modelForMatrices$Z)
    Zu <- vector(mode = "list", length = nz) # list for Zu
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
      # print(originalModel$sigma)
      G.mv.original.List <- lapply(as.list((zToUse)),function(x){
        kronecker(as.matrix(originalModelForMatricesSE$K[[x]]),originalModel$sigma[,,x])
      }); G.mv.original <- do.call(adiag1,G.mv.original.List)
      
      tZ.mv.original.List <- lapply(as.list((zToUse)),function(x){
        kronecker(as.matrix(t(originalModelForMatricesSE$Z[[x]])),TT)
      }); tZ.mv.original <- do.call(rbind,tZ.mv.original.List)
      
      G.tZ.mv.original <- G.mv.original %*% tZ.mv.original # GZ' as many rows as obs, as many cols as levels
      
      if(length(fToUse) > 0){ # Cov(b,u) = 0 - (X'ViX)-X' Vi GZ
        cov.b.pev <- 0 - ( XtViX %*% t(X.mv.original) %*% originalModel$Vi %*% t(G.tZ.mv.original) )
        GZtViZG <-  G.tZ.mv.original%*%originalModel$Vi%*%t(G.tZ.mv.original)
        GZtViXPXViZG <-  G.tZ.mv.original %*% originalModel$Vi %*% X.mv.original %*% XtViX %*% t(X.mv.original) %*% originalModel$Vi %*% t(G.tZ.mv.original)
        VarU <- GZtViZG - GZtViXPXViZG
        pev <- G.mv.original - VarU
      }
        
      # (185 x 3)' (185 x 185) = (3 x 185)  (185 x 185) (levs164 x obs185)' = (3 x 164)
      # bring the design matrices for extended model to get standard errors for the extended model
      Z.extended <- lapply(as.list((zToUse)),function(x){
        kronecker(as.matrix(t(modelForMatrices$Z[[x]])),TT)
      }); Z.extended <- do.call(rbind,Z.extended)
      tZ.extended <- t(Z.extended)
    }
  }
  ##################################################
  # add them up
  # print(zToUse)
  if(length(fToUse) > 0 & !is.null(zToUse)){ # fixed and random effects included
    y.hat <- Xb + Reduce("+",Zu[zToUse]) # y.hat = Xb + Zu.1 + ... + Zu.n
    standard.errors <- sqrt(abs(
                              rowSums((X.mv.extended%*%XtViX)*X.mv.extended) +
                              rowSums(2*(X.mv.extended%*%cov.b.pev)*tZ.extended) +
                              rowSums((tZ.extended%*%pev)*tZ.extended)
                               )
                            )
  }else if( length(fToUse) == 0 & !is.null(zToUse) ){ # only random effects included
    y.hat <- Reduce("+",Zu[zToUse]) # y.hat = Xb + Zu.1 + ... + Zu.n
    standard.errors <- sqrt(abs(
                                rowSums((tZ.extended%*%pev)*tZ.extended)
                                )
                            )
    XtViX=NULL
  }else if( length(fToUse) > 0 & is.null(zToUse) ){ # only fixed effects included
    y.hat <- Xb # y.hat = Xb 
    standard.errors <- sqrt(abs(
                                rowSums((X.mv.extended%*%XtViX)*X.mv.extended)
                                )
                            )
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
  
  myForm <- paste0("predicted.value~",paste(classify, collapse = "+"), "+trait")
  pvals <- aggregate(as.formula(myForm), FUN=mean, data=DTX2)
  myFormSE <- paste0("standard.error~",paste(classify, collapse = "+"), "+trait")
  pvalsSE <- aggregate(as.formula(myFormSE), FUN=mean, data=DTX2)
  # pvalsSE <- aggregate(as.formula(myForm), FUN=function(x){1.96 * (sd(x)/sqrt(length(x)))}, data=DTX2)
  colnames(pvalsSE)[ncol(pvalsSE)] <- "standard.error"
  # print(classify)
  pvals <- merge(pvals,pvalsSE, by=c("trait",classify))
  
  # ##################################################
  ## hypertable summary
  nLevels <- c(unlist(lapply(originalModelForMatricesSE$X,ncol)), unlist(lapply(originalModelForMatricesSE$Z,ncol)))
  namesLevels <- c(unlist(oto$fixed),names(object$U))
  
  namesLevelsO <- data.frame( x=c(unlist(oto2$fixed),unlist(oto2$random)), y=c(object$termsN$fixed, object$termsN$random))
  namesLevelsO <- unlist(apply(namesLevelsO,1,function(x){rep(x[1],x[2])}))
  
  formLevels <- c(rep("fixed",length(unlist(oto$fixed))),rep("random",length(names(object$U))))
  id <- 1:length(formLevels)
  ignored <- rep(TRUE,length(id)) # all ignored by default
  include <- rep(FALSE,length(id)) # none include by default
  
  include[fToUse]=TRUE # specify which fixed effects where used
  if(!is.null(zToUse)){include[length(modelForMatrices$X)+zToUse]=TRUE} # specify which random effects where used
  
  # print(fToUse)
  # print(length(modelForMatrices$X)+zToUse)
  ignored[fToUse]=FALSE # specify which fixed effects ARE NOT ignored
  if(!is.null(zToUse)){ignored[length(modelForMatrices$X)+zToUse]=FALSE} # specify which random effects ARE NOT ignored
  
  predictSummary <- data.frame(namesLevelsO,namesLevels,formLevels,nLevels,id,ignored,include)
  colnames(predictSummary) <- c("termHL","term","type","nLevels","id","ignored","include")
  
  toreturn2 <- list(pvals=pvals,
                    hypertable=predictSummary,
                    model=modelForMatrices,
                    C11=XtViX,
                    C12=cov.b.pev,
                    C22=pev
                    
  )
  attr(toreturn2, "class")<-c("predict.mmer", "list")
  
  return(toreturn2)
}


"print.predict.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("
    The predictions are obtained by averaging/aggregating across 
    the hypertable calculated from model terms constructed solely 
    from factors in the include sets. You can customize the model 
    terms used with the 'hypertable' argument. Current model terms used:\n")
  ))
  print(x$hypertable)
  cat(blue(paste("\n Head of predictions:\n")
  ))
  head(x$pvals,...)
}

"head.predict.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste(" 
    The predictions are obtained by averaging/aggregating across 
    the hypertable calculated from model terms constructed solely 
    from factors in the include sets. You can customize the model 
    terms used with the 'hypertable' argument. Current model terms used:\n")
  ))
  print(x$hypertable)
  cat(blue(paste("\n Head of predictions:\n")
  ))
  head(x$pvals,...)
}


"tail.predict.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("
    The predictions are obtained by averaging/aggregating across 
    the hypertable calculated from model terms constructed solely 
    from factors in the include sets. You can customize the model 
    terms used with the 'hypertable' argument. Current model terms used:\n")
  ))
  print(x$hypertable)
  cat(blue(paste("\n Tail of predictions:\n")
  ))
  tail(x$pvals,...)
}

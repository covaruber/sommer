
#### =========== ######
## PREDICT FUNCTION #
#### =========== ######
"predict.mmer" <- function(object,classify=NULL,RtermsToForce=NULL,FtermsToForce=NULL,...){
  
  newdata=NULL
  if(is.null(newdata)){
    newdata <- object$data
  }
  
  if(is.null(object$call$random)){
    prov <- mmer(fixed=object$call$fixed,
                 #random=object$call$random,
                 rcov=object$call$rcov,
                 data=object$data, return.param = TRUE,#reshape.results=TRUE,
                 na.method.Y = object$call$na.method.Y, 
                 na.method.X = object$call$na.method.X)
    prov2 <- mmer(fixed=object$call$fixed,
                  # random=object$call$random,
                  rcov=object$call$rcov, iters=1,
                  data=newdata, return.param = FALSE,reshape.output=FALSE,
                  init = object$sigma_scaled, constraints = object$constraints,
                  na.method.Y = object$call$na.method.Y, 
                  na.method.X = object$call$na.method.X)
  }else{
    prov <- mmer(fixed=object$call$fixed,
                 random=object$call$random,
                 rcov=object$call$rcov,
                 data=object$data, return.param = TRUE,#reshape.results=TRUE,
                 na.method.Y = object$call$na.method.Y, 
                 na.method.X = object$call$na.method.X)
    prov2 <- mmer(fixed=object$call$fixed,
                  random=object$call$random,
                  rcov=object$call$rcov,
                  data=newdata, return.param = FALSE,reshape.output =FALSE,
                  init = object$sigma_scaled, constraints = object$constraints,
                  na.method.Y = object$call$na.method.Y, iters = 1,
                  na.method.X = object$call$na.method.X)
  }
  
  nt <- ncol(prov[[1]])
  
  if(is.null(classify)){
    # cat(paste("Returning predictions including all random effects.\n"))
    classify.split <- names(object$U)
  }else{
    classify.split <- unique(c(classify,unlist(strsplit(classify,":"))))
    classify.split <- unique(c(classify,classify.split))
    # classify.split <- classify
  }
  ## term names in the random effects
  if(is.null(classify) & is.null(RtermsToForce)){
    RtermsToForce <- 1:length(object$U)
  }
  # if(include=="all"){
  if(!is.null(object$call$random)){
    namesres <- lapply(as.list(names(object$U)),function(x){
      strsplit(x,":")[[1]]
    })
    ## the != 0 should be used for the predict
    presence <- unlist(lapply(namesres,function(x){
      length(which(x %in% classify.split))
    }))
    arepresent <- which(presence!=0)
    if(!is.null(RtermsToForce)){
      arepresent <- RtermsToForce
    }
    used <- unlist(names(object$U))[arepresent]
    # cat(paste("Random effects included in the predictions:\n",used))
    # cat("\n")
  }else{ #there's no random effects
    arepresent <- integer()
    used <- character()
  }
  ## get the index for each fixed effect
  nfixedeff <- length(prov[[2]])
  if(nfixedeff != length(prov[[19]])){
    fnames <- c("(Intercept)",prov[[19]])
  }else{fnames <- c(prov[[19]])}
  
  fnamesindex <- list()
  counter <- 1
  kept <- numeric()
  for(i in 1:nfixedeff){
    if(ncol(prov[[2]][[i]]) > 0){ # the fixed effect was fitted
      pX <- kronecker(prov[[2]][[i]],prov[[3]][[i]])
      # pX <- kronecker(prov[[3]][[i]],prov[[2]][[i]])
      endcounter <- (counter+ncol(pX)-1)
      fnamesindex[[i]] <- counter:endcounter
      counter <- endcounter+1
      kept[i] <- i
    }
  }
  fnamesindex <- fnamesindex[which(unlist(lapply(fnamesindex,length)) > 0)]
  nfixedeff <- length(fnamesindex)
  fnames <- fnames[na.omit(kept)]
  # fpicked <- which(fnames %in% c("(Intercept)",classify.split))
  if(!is.null(FtermsToForce)){
    fpicked <- intersect(FtermsToForce,1:nfixedeff)
  }else{
    fpicked <- 1:length(fnames)
  }
  usedf <- fnames[fpicked] # used fixed
  # cat(paste("Fixed effects included in the predictions:\n",usedf))
  # cat("\n")
  fpickedindex <- fnamesindex[fpicked]
  fpickedindex <- unlist(fpickedindex)
  
  fnonpicked <- setdiff(1:length(fnamesindex),fpicked)
  
  nonusedf <- fnames[fnonpicked] # simple averaging set
  # fnonpicked <- setdiff(1:length(fnames),fpicked)
  # fnonpickedindex <- fnamesindex[fnonpicked]
  # fnonpickedindex <- unlist(fnonpickedindex)
  # print(fpickedindex)
  # print(fnonpickedindex)
  ## build the matrices
  if(length(arepresent) > 0){
    ## build the multivariate K and Z
    VarKl <- lapply(as.list((arepresent)),function(x){
      kronecker(as.matrix(prov[[5]][[x]]),object$sigma[[x]])
    })
    VarK <- do.call(adiag1,VarKl)
    
    Zl <- lapply(as.list((arepresent)),function(x){
      kronecker(as.matrix(t(prov[[4]][[x]])),diag(nt))
    })
    tZm <- do.call(rbind,Zl)
    
    ZKfv <- VarK %*% tZm
    
    Xlist <- list()
    for(o in 1:length(prov[[3]])){
      Xlist[[o]] <- kronecker(prov[[2]][[o]],prov[[3]][[o]])
    }
    Xm <- do.call(cbind,Xlist)
    Xmv <- Xm
    if(length(fnonpicked) > 0){ #then we have to average across those terms
      for(ifnp in fnonpicked){
        fnonpickedindex <- fnamesindex[[ifnp]]
        if(ncol(as.matrix(Xm[,fnonpickedindex])) > 1){
          Xm[,fnonpickedindex] <- ((Xm[,fnonpickedindex]*0)+1)/(ncol(as.matrix(Xm[,fnonpickedindex]))+1)
        }else{
          Xm[,fnonpickedindex] <- mean(Xm[,fnonpickedindex],na.rm=TRUE)
        }
      }
    }
    Xmv <- Xm# Xmv/ncol(Xmv)
    # Xm <- Xm/ncol(Xm)
    
    #### cov(b,u - u.hat), xvxi, pev (C12,C11,C22)
    cov.b.pev <- 0 - t(prov2$P %*% (Xmv)) %*% prov2$Vi %*% t(ZKfv)
    # xvxi <- as.matrix(prov2$VarBeta[fpickedindex,fpickedindex])
    xvxi <- as.matrix(prov2$VarBeta)
    # print(dim(Xm))
    # print(dim(xvxi))
    # xvxi <- xvxi*3
    # xvxi[2:3,2:3] <- xvxi[2:3,2:3]/3
    pev <- do.call(adiag1,prov2$PevU[arepresent])
    
    MMsp <- Xm# for predictions
    MMspv <- Xmv# for SE's
    Zl2 <- lapply(as.list((arepresent)),function(x){
      kronecker(as.matrix(t(prov[[4]][[x]])),diag(nt))
    })
    
    tZm2 <- do.call(rbind,Zl2)
    MMnsp <- t(tZm2) ## random term
    
    standard.errors <- sqrt(rowSums((MMspv%*%xvxi)*MMspv) + 
                              rowSums(2*(MMspv%*%cov.b.pev)*MMnsp) + 
                              rowSums((MMnsp%*%pev)*MMnsp))
    
    blups <- unlist(prov2$U[arepresent])
    blues <- prov2$Beta#[fpickedindex,]
    
    coeffs <- c(blues, blups)
    # 
    predicted.value <- as.vector(cbind(MMsp, MMnsp)%*%as.vector(coeffs))
    
    newd <- as.data.frame(newdata)#data.frame(,predicted.value,standard.errors)
    
    preds <- as.data.frame(matrix(predicted.value,ncol=nt,byrow = TRUE))
    ses <- as.data.frame(matrix(standard.errors,ncol=nt,byrow = TRUE))
    colnames(preds) <- paste("predicted.value",colnames(prov[[1]]),sep=".")
    colnames(ses) <- paste("standard.errors",colnames(prov[[1]]),sep=".")
    
    newd <- data.frame(newd,preds,ses)
    
  }else{
    
    Xlist <- list()
    for(o in 1:length(prov[[3]])){
      Xlist[[o]] <- kronecker(prov[[2]][[o]],prov[[3]][[o]])
    }
    Xm <- do.call(cbind,Xlist)
    # Xmv <- Xm
    # print(fnonpicked)
    if(length(fnonpicked) > 0){ #then we have to average across those terms
      for(ifnp in fnonpicked){ # ifnp <- fnonpicked[5]
        fnonpickedindex <- fnamesindex[[ifnp]]
        if(ncol(as.matrix(Xm[,fnonpickedindex])) > 1){ ## if is factor or character type
          Xm[,fnonpickedindex] <- ((Xm[,fnonpickedindex]*0)+1)/(ncol(as.matrix(Xm[,fnonpickedindex]))+1)
        }else{ # if is a numeric covariate we average across the covariate
          Xm[,fnonpickedindex] <- mean(Xm[,fnonpickedindex],na.rm=TRUE)
        }
      }
    }
    # print(head(Xm))
    Xmv <- Xm#Xmv/ncol(Xmv)
    
    xvxi <- as.matrix(prov2$VarBeta)
    MMsp <- Xm## fixed part
    MMspv <- Xmv## fixed part
    standard.errors <- sqrt(rowSums((MMspv%*%xvxi)*MMspv))
    coeffs <- c(prov2$Beta) #[fpickedindex,]
    predicted.value <- as.vector(cbind(MMsp)%*%as.vector(coeffs))
    newd <- as.data.frame(newdata)
    preds <- as.data.frame(matrix(predicted.value,ncol=nt,byrow = TRUE))
    ses <- as.data.frame(matrix(standard.errors,ncol=nt,byrow = TRUE))
    colnames(preds) <- paste("predicted.value",colnames(prov[[1]]),sep=".")
    colnames(ses) <- paste("standard.errors",colnames(prov[[1]]),sep=".")
    newd <- data.frame(newd,preds,ses)
    cov.b.pev <- NULL
    pev <- NULL
  }
  
  if(is.null(classify)){
    toreturn2 <- list(predictions=newd,RE.used=used, FE.used=usedf,fitted=newd, 
                      C11=xvxi, C12=cov.b.pev, C22=pev)
  }else{
    resp0 <- paste(c(colnames(preds),colnames(ses)),collapse = ",")
    myform <- paste(paste("cbind(",resp0,")~"),paste(classify.split,collapse = "+"))
    toreturn <- aggregate(as.formula(myform),data=newd, FUN=mean,na.rm=TRUE)
    toreturn2 <- list(predictions=toreturn,RE.used=used, FE.used=usedf,
                      FE.nonused=nonusedf,fitted=newd, 
                      C11=xvxi, C12=cov.b.pev, C22=pev)
  }
  
  attr(toreturn2, "class")<-c("predict.mmer", "list")
  
  return(toreturn2)
}

"print.predict.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  aver <- paste(x$FE.nonused,collapse = ", ")
  fused <- paste(x$FE.used,collapse = ", ")
  rused <- paste(x$RE.used,collapse = ", ")
  cat(blue(paste("\n  The predictions are obtained by averaging across the hypertable
  calculated from model terms constructed solely from factors in
  the averaging and classify sets.
 - Use 'FtermsToForce' to move ignored fixed factors into the averaging set.
 - Use 'RtermsToForce' to move ignored random factors into the averaging set.
 - The simple averaging set:",aver,"\n",
  "- The fixed effects included:",fused,"\n",
  "- The random effects included:",rused,"\n\n")
  ))
  head(x$predictions,...)
  # head(print.predict.mmer(pp))
  # print((x$predictions))
}

"head.predict.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  aver <- paste(x$FE.nonused,collapse = ", ")
  fused <- paste(x$FE.used,collapse = ", ")
  rused <- paste(x$RE.used,collapse = ", ")
  cat(blue(paste("\n  The predictions are obtained by averaging across the hypertable
  calculated from model terms constructed solely from factors in
  the averaging and classify sets.
 - Use 'FtermsToForce' to move ignored fixed factors into the averaging set.
 - Use 'RtermsToForce' to move ignored random factors into the averaging set.
 - The simple averaging set:",aver,"\n",
  "- The fixed effects included:",fused,"\n",
  "- The random effects included:",rused,"\n\n")
  ))
  head(x$predictions, ...)
  # head(print.predict.mmer(pp))
  # print((x$predictions))
}

#### =========== ####
## SUMMARY FUNCTION mmer #
#### =========== ####
"summary.mmer" <- function(object, ...) {
  
  #dim(object$u.hat)
  digits = max(3, getOption("digits") - 3)
  #forget <- length(object)
  
  groupss.nn <- lapply(object$U,function(x){
    unlist(lapply(x,length))
  })
  groupss.nn <- do.call(rbind,groupss.nn)
  
  lll <- object$monitor[1,]
  lll <- lll[which(lll > 1e-300 | lll < 0)]
  lll2 <- lll[length(lll)]
  
  LLAIC <- data.frame(as.numeric(lll2), as.numeric(object$AIC),
                      as.numeric(object$BIC), object$method, object$convergence)
  colnames(LLAIC) = c("logLik","AIC","BIC","Method","Converge")
  rownames(LLAIC) <- "Value"
  
  method=object$method
  #extract fixed effects
  coef <- as.data.frame((object$Beta))#, Std.Error=(matrix(sqrt(diag(object$Var.beta.hat)),ncol=1)), t.value=(matrix((object$beta.hat-0)/sqrt(diag(object$Var.beta.hat)), ncol=1)))
  # if(dim(coef)[1] == 1){rownames(coef) <- "Intercept"}
  
  ## se and t values for fixed effects
  ts <- ncol(object$sigma[[1]])
  s2.beta <- diag(as.matrix(object$VarBeta))
  coef$Std.Error <- sqrt(abs(s2.beta))
  coef$t.value <- coef$Estimate/coef$Std.Error
  # print(coef)
  # nse.beta <- length(s2.beta)/ts
  # inits <- seq(1,length(s2.beta),nse.beta)
  # ends <- inits+nse.beta-1
  # seti <- list() # stardard errors partitioned by trait
  # for(u in 1:ts){
  #   prox <- data.frame(coef[,u],sqrt(abs(s2.beta[inits[u]:ends[u]])))
  #   prox$`t value` <- prox[,1]/prox[,2]
  #   colnames(prox) <- c("Estimate","Std. Error","t value")
  #   rownames(prox) <- rownames(coef)
  #   seti[[u]] <- prox
  # }
  # names(seti) <- colnames(object$sigma[[1]])
  
  vcsl <- list()
  consl <- list()
  for(k in 1:length(object$sigma)){
    x <- object$sigma[[k]]
    y <- object$constraints[[k]]
    xn <- names(object$sigma)[k]
    vcs <- numeric()
    cons <- numeric()
    counter <-1
    for(i in 1:ncol(x)){
      for(j in i:ncol(x)){
        # print(y[i,j])
        if( y[i,j] != 0 ){
          vcs[counter] <- x[i,j]
          cons[counter] <- y[i,j]
          names(vcs)[counter] <- paste(colnames(x)[i],colnames(x)[j],sep="-" )
          counter <- counter+1
        }
      }
    }
    vcsl[[xn]] <- as.data.frame(vcs)
    consl[[xn]] <- as.data.frame(cons)
  }
  mys2 <- do.call(rbind,vcsl)
  mycons <- do.call(rbind,consl)
  
  rrr <- lapply(vcsl,rownames)
  rrr <- rrr[which(unlist(lapply(rrr, length)) > 0)]
  for(o in 1:length(rrr)){rrr[[o]] <- paste(names(rrr)[o],rrr[[o]],sep=".")}
  rownames(mys2) <- as.vector(unlist(rrr))
  
  varcomp <- as.data.frame(cbind(mys2,sqrt(abs(diag(object$sigmaSE)))))
  varcomp[,3] <- varcomp[,1]/varcomp[,2]
  colnames(varcomp) <- c("VarComp","VarCompSE","Zratio")
  varcomp$Constraint <- replace.values(mycons$cons, 1:3, c("Positive","Unconstr","Fixed"))
  
  output <- list(groups=groupss.nn, varcomp=varcomp, betas=coef, method=method,logo=LLAIC)
  attr(output, "class")<-c("summary.mmer", "list")
  return(output)
}

"print.summary.mmer"<-function (x, digits = max(3, getOption("digits") - 3),  ...){
  
  nmaxchar0 <- max(as.vector(unlist(apply(data.frame(rownames(x$varcomp)),1,nchar))),na.rm = TRUE)
  
  if(nmaxchar0 < 26){
    nmaxchar0 <- 26
  } # + 26 spaces we have nmaxchar0+26  spaces to put the title
  
  nmaxchar <- nmaxchar0+34 ## add spaces from the 3 columns
  nmaxchar2 <- nmaxchar0+18
  nmaxchar3 <- nmaxchar0+34-46 #round(nmaxchar0/2)
  rlh <- paste(rep("*",round(nmaxchar2/2)),collapse = "")
  rlt <- paste(rep(" ",ceiling(nmaxchar3/2)),collapse = "")
  digits = max(3, getOption("digits") - 3)
  ################################################
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat(paste("\n",rlt,"Multivariate Linear Mixed Model fit by REML",rlt,"\n", collapse = ""))
  cat(paste(rlh," sommer 4.0 ",rlh, "\n", collapse = ""))
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\n")
  cat("")
  print(x$logo)#, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nVariance-Covariance components:\n")
  print(x$varcomp, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nFixed effects:\n")
  if(nrow(x$betas) > 8){
    print(x$betas[1:8,], digits = digits)
    cat("   ... please access the object to see more\n")
  }else{
    print(x$betas, digits = digits)
  }
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nGroups and observations:\n")
  print(x$groups, digits = digits)
  cat(paste(rep("=",nmaxchar), collapse = ""))
  cat("\nUse the '$' sign to access results and parameters")#\nArguments set to FALSE for multiresponse models:\n'draw', and 'gwas.plots'\n")
  ################################################
}

#### =========== ####
## SUMMARY FUNCTION mmergwas #
#### =========== ####
# "summary.mmergwas" <- function(object, ...) {
#   
#   #dim(object$u.hat)
#   digits = max(3, getOption("digits") - 3)
#   #forget <- length(object)
#   
#   groupss.nn <- lapply(object$U,function(x){
#     dim(x)
#   })
#   groupss.nn <- as.data.frame(do.call(rbind,groupss.nn))
#   
#   lll <- object$monitor[1,]
#   lll <- lll[which(lll > 1e-300 | lll < 0)]
#   lll2 <- lll[length(lll)]
#   
#   LLAIC <- data.frame(as.numeric(lll2), as.numeric(object$AIC),
#                       as.numeric(object$BIC), object$method, object$convergence)
#   colnames(LLAIC) = c("logLik","AIC","BIC","Method","Converge")
#   rownames(LLAIC) <- "Value"
#   
#   method=object$method
#   #extract fixed effects
#   coef <- data.frame(Estimate=unlist(object$Beta))#, Std.Error=(matrix(sqrt(diag(object$Var.beta.hat)),ncol=1)), t.value=(matrix((object$beta.hat-0)/sqrt(diag(object$Var.beta.hat)), ncol=1)))
#   # if(dim(coef)[1] == 1){rownames(coef) <- "Intercept"}
#   
#   ## se and t values for fixed effects
#   ts <- dim(object$sigma)[1]
#   s2.beta <- diag(as.matrix(object$VarBeta))
#   coef$Std.Error <- sqrt(abs(s2.beta))
#   coef$t.value <- coef$Estimate/coef$Std.Error
#   # print(coef)
#   # nse.beta <- length(s2.beta)/ts
#   # inits <- seq(1,length(s2.beta),nse.beta)
#   # ends <- inits+nse.beta-1
#   # seti <- list() # stardard errors partitioned by trait
#   # for(u in 1:ts){
#   #   prox <- data.frame(coef[,u],sqrt(abs(s2.beta[inits[u]:ends[u]])))
#   #   prox$`t value` <- prox[,1]/prox[,2]
#   #   colnames(prox) <- c("Estimate","Std. Error","t value")
#   #   rownames(prox) <- rownames(coef)
#   #   seti[[u]] <- prox
#   # }
#   # names(seti) <- colnames(object$sigma[[1]])
#   
#   vcsl <- list()
#   consl <- list()
#   nnn <- dim(object$sigma)
#   for(k in 1:nnn[3]){
#     x <- as.matrix(object$sigma[,,k])
#     y <- object$constraints[[k]]
#     xn <- k#names(object$sigma)[k]
#     vcs <- numeric()
#     cons <- numeric()
#     counter <-1
#     for(i in 1:ncol(x)){
#       for(j in i:ncol(x)){
#         # print(y[i,j])
#         if( y[i,j] != 0 ){
#           vcs[counter] <- x[i,j]
#           cons[counter] <- y[i,j]
#           names(vcs)[counter] <- paste(colnames(x)[i],colnames(x)[j],sep="-" )
#           counter <- counter+1
#         }
#       }
#     }
#     vcsl[[xn]] <- as.data.frame(vcs)
#     consl[[xn]] <- as.data.frame(cons)
#   }
#   mys2 <- do.call(rbind,vcsl)
#   mycons <- do.call(rbind,consl)
#   
#   # rrr <- lapply(vcsl,rownames)
#   # rrr <- rrr[which(unlist(lapply(rrr, length)) > 0)]
#   # for(o in 1:length(rrr)){rrr[[o]] <- paste(names(rrr)[o],rrr[[o]],sep=".")}
#   # rownames(mys2) <- as.vector(unlist(rrr))
#   
#   varcomp <- as.data.frame(cbind(mys2,sqrt(abs(diag(object$sigmaSE)))))
#   varcomp[,3] <- varcomp[,1]/varcomp[,2]
#   colnames(varcomp) <- c("VarComp","VarCompSE","Zratio")
#   varcomp$Constraint <- replace.values(mycons$cons, 1:3, c("Positive","Unconstr","Fixed"))
#   
#   output <- list(groups=groupss.nn, varcomp=varcomp, betas=coef, method=method,logo=LLAIC)
#   attr(output, "class")<-c("summary.mmergwas", "list")
#   return(output)
# }
# 
# "print.summary.mmergwas"<-function (x, digits = max(3, getOption("digits") - 3),  ...){
#   
#   nmaxchar0 <- max(as.vector(unlist(apply(data.frame(rownames(x$varcomp)),1,nchar))),na.rm = TRUE)
#   
#   if(nmaxchar0 < 26){
#     nmaxchar0 <- 26
#   } # + 26 spaces we have nmaxchar0+26  spaces to put the title
#   
#   nmaxchar <- nmaxchar0+34 ## add spaces from the 3 columns
#   nmaxchar2 <- nmaxchar0+18
#   nmaxchar3 <- nmaxchar0+34-46 #round(nmaxchar0/2)
#   rlh <- paste(rep("*",round(nmaxchar2/2)),collapse = "")
#   rlt <- paste(rep(" ",ceiling(nmaxchar3/2)),collapse = "")
#   digits = max(3, getOption("digits") - 3)
#   ################################################
#   cat(paste(rep("=",nmaxchar), collapse = ""))
#   cat(paste("\n",rlt,"Multivariate Linear Mixed Model fit by REML",rlt,"\n", collapse = ""))
#   cat(paste(rlh," sommer 4.0 ",rlh, "\n", collapse = ""))
#   cat(paste(rep("=",nmaxchar), collapse = ""))
#   cat("\n")
#   cat("")
#   print(x$logo)#, digits = digits)
#   cat(paste(rep("=",nmaxchar), collapse = ""))
#   cat("\nScaled Variance-Covariance components:\n")
#   print(x$varcomp, digits = digits)
#   cat(paste(rep("=",nmaxchar), collapse = ""))
#   cat("\nScaled Fixed effects:\n")
#   if(nrow(x$betas) > 8){
#     print(x$betas[1:8,], digits = digits)
#     cat("   ... please access the object to see more\n")
#   }else{
#     print(x$betas, digits = digits)
#   }
#   cat(paste(rep("=",nmaxchar), collapse = ""))
#   cat("\nGroups and observations:\n")
#   print(x$groups, digits = digits)
#   cat(paste(rep("=",nmaxchar), collapse = ""))
#   cat("\nUse the '$' sign to access results and parameters")#\nArguments set to FALSE for multiresponse models:\n'draw', and 'gwas.plots'\n")
#   ################################################
# }

#### =========== ######
## RESIDUALS FUNCTION #
#### =========== ######
"residuals.mmer" <- function(object, ...) {
  digits = max(3, getOption("digits") - 3)
  
  # if(type=="conditional"){
  #   output <- object$cond.residuals
  #   #colnames(output) <- names(object)
  # }else{
  output <- object$residuals
  #colnames(output) <- names(object)
  # }
  return(output)
}

"print.residuals.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
}
#### =========== ######
## RANEF FUNCTION #
#### =========== ######

"randef" <- function(object) {
  # if(class(object)=="mmer"){
  #   digits = max(3, getOption("digits") - 3)
  #   cat("Returning object of class 'list' where each element correspond to one random effect.")
  #   output <- object$u.hat
  # }else{
  #   stop("Class not recognized.\n",call. = FALSE)
  # }
  output<- object$U
  return(output)
}

#"print.ranef.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
#  print((x))
#}

#### =========== ######
## FIXEF FUNCTION #
#### =========== ######


#"fixef.mmer" <- function(object, ...) {
#  digits = max(3, getOption("digits") - 3)
#  output <- object$beta.hat
#  return(output)
#}

#"print.fixef.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
#  print((x))
#}

#### =========== ####
## FITTED FUNCTION ##
#### =========== ####

"fitted.mmer" <- function(object, type="complete", ...) {
  #type="complete" 
  # digits = max(3, getOption("digits") - 3)
  # if(type=="complete"){
  #   output<- object$fitted.y
  #   #colnames(output) <- names(object)
  # }else{
  #   output<- object$fitted.u
  #   #colnames(output) <- names(object)
  # }
  output<- object$fitted
  return(output)
}

"print.fitted.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
} 


#### =========== ####
## COEF FUNCTION ####
#### =========== ####

"coef.mmer" <- function(object, ...){
  object$Beta
}

"print.coef.mmer"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print((x))
} 

#### =========== ####
## ANOVA FUNCTION ###
#### =========== ####
anova.mmer <- function(object, object2=NULL, type=1, ...) {
  sequential.fit <- function(object, type=1){
    
    fixed <- object$call$fixed
    yuyuf <- strsplit(as.character(fixed[3]), split = "[+-]")[[1]]
    fixedtermss <- apply(data.frame(yuyuf),1,function(x){
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    fixedtermss <- fixedtermss[which(fixedtermss != "-1")]
    fixedtermss <- fixedtermss[which(fixedtermss != "1")]
    fixedtermss <- c("1",fixedtermss)
    response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
    # responsef <- as.formula(paste(response,"~1"))
    # fixedtermss <- as.list(fixedtermss)
    # fixedtermss[[length(fixedtermss)+1]] <- unlist(fixedtermss)
    SSM <- numeric(length = length(fixedtermss)+1)
    df <- numeric(length = length(fixedtermss)+1)
    models <- character(length = length(fixedtermss)+1)
    # if(length(fixedtermss) < 0){
    #   until <- 1
    # }else{
    until <- length(fixedtermss)
    # }
    
    namesl <- character(length = length(fixedtermss)+1)
    for(i in (until+1):1){
      
      if(i == (length(fixedtermss)+1)){ # first the full model
        myf <- paste(response,"~",paste(fixedtermss,collapse = " + "))
        fixedi <- as.formula(myf)
        usef <- fixedtermss
      }else if(i == 1){
        myf <- paste(response,"~1")
        fixedi <- as.formula(myf)
      }else{
        if(type == 1){
          usef <- setdiff(fixedtermss,fixedtermss[i])
          myf <- paste(response,"~",paste(usef,collapse = " + "))
          fixedi <- as.formula(myf)
        }else if(type == 2){
          usef <- fixedtermss[1:i]
          myf <- paste(response,"~",paste(usef,collapse = " + "))
          fixedi <- as.formula(myf)
        }else{
          stop("anova type not enabled.", call. = FALSE)
        }
      }
      namesl[i] <- paste(usef,collapse = ",")
      models[i] <- myf
      cat("\n")
      cat(paste("Fitting:", myf))
      # cat("\n")
      # print(fixedi)
      if(is.null(object$call$random)){
        prov0 <- mmer(fixed=fixedi,
                      # random = object$call$random,
                      rcov=object$call$rcov, iters=1,
                      data=object$data, return.param = FALSE,reshape.output=FALSE,
                      verbose = FALSE,
                      init = object$sigma_scaled, constraints = object$constraints,
                      na.method.Y = object$call$na.method.Y, 
                      na.method.X = object$call$na.method.X)
        prov <- mmer(fixed=fixedi,
                     # random = object$call$random,
                     rcov=object$call$rcov,
                     data=object$data, return.param = TRUE,
                     na.method.Y = object$call$na.method.Y, 
                     na.method.X = object$call$na.method.X)
      }else{
      prov0 <- mmer(fixed=fixedi,
                    random = object$call$random,
                    rcov=object$call$rcov, iters=1,
                    data=object$data, return.param = FALSE,reshape.output=FALSE,
                    verbose = FALSE,
                    init = object$sigma_scaled, constraints = object$constraints,
                    na.method.Y = object$call$na.method.Y, 
                    na.method.X = object$call$na.method.X)
      prov <- mmer(fixed=fixedi,
                   random = object$call$random,
                   rcov=object$call$rcov,
                   data=object$data, return.param = TRUE,
                   na.method.Y = object$call$na.method.Y, 
                   na.method.X = object$call$na.method.X)
      }
      Y <- prov[[1]]
      n= nrow(Y)
      J = matrix(1, nrow=n, ncol=n)
      # Total sum of Square 
      SSTO = t(Y) %*% Y - (1/n)*t(Y)%*%J%*%Y
      Xlist <- list()
      # ncolsx <- numeric()
      for(o in 1:length(prov[[3]])){
        Xlist[[o]] <- kronecker(prov[[2]][[o]],prov[[3]][[o]])
        # ncolsx[o] <- ncol(Xlist[[o]])
      }
      X <- do.call(cbind,Xlist)
      head(X)
      ncolsx <- unlist(lapply(Xlist,ncol))
      ncolsx <- ncolsx[which(ncolsx>0)]
      xdf <- ncolsx[length(ncolsx)]
      #X <- as.matrix(X[,1:2])
      b = as.matrix(prov0$Beta)
      # regression sum of square 
      ss <- (SSR = t(b)%*%t(X)%*%Y - (1/n)%*%t(Y)%*%J%*%Y)
      if(i == (length(fixedtermss)+1)){
        SSM[i] = as.numeric(ss)
      }else{
        SSM[i] = SSM[length(fixedtermss)+1] - as.numeric(ss)
      }
      df[i] <- xdf#ncol(X)
    }
    
    myanova <- data.frame(Df=rev(df), Sum.Sq=rev(SSM))
    
    # print(myanova)
    # print(c("Full",rev(namesl[-c(length(namesl))])))
    # rownames(myanova) <- c("Full",rev(namesl[-c(length(namesl))]))
    
   #else{
    #   rownames(myanova) <- c("Full", "None", rev(fixedtermss[-c(1)]))
    # }
    myanova$Mean.Sq <- myanova$Sum.Sq/myanova$Df
    
    vc <- object$sigma
    vare <- as.numeric(vc[[length(vc)]])
    dfe <- nrow(Y)-sum(myanova$Df)
    myanova2 <- data.frame(Df=dfe, Sum.Sq=NA,Mean.Sq=vare)
    rownames(myanova2) <- "Residuals"
    myanovaf <- rbind(myanova,myanova2)
    
    myanovaf$F.value <- myanovaf$Mean.Sq/vare
    # print(c(myanovaf,dfe))
    myanovaf$`Pr(>F)` <- round(apply(myanovaf,1,function(x){1-pf(x[4],x[1],dfe)}),4)
    myanovaf[nrow(myanovaf),4:5] <- NA
    myanovaf$Models <- c(rev(models),"")
    if(type==1){
      rownames(myanovaf) <- c("Full", rev(fixedtermss), "Residuals")
      myanovaf <- myanovaf[-c(nrow(myanovaf)-1),]
    }else{
      rownames(myanovaf) <- c("Full","all", rev(fixedtermss[-c(1)]), "Residuals")
      # print(rev(fixedtermss))
    }
    cat("\n\n")
    return(myanovaf)
  }
  signifo <- function(x){
    if(x >= 0 & x < 0.001){y="***"}
    if(x >= 0.001 & x < 0.01){y="**"}
    if(x >= 0.01 & x < 0.05){y="*"}
    if(x >= 0.05 & x < 0.1){y="."}
    if(x > 0.1){y=""}
    return(y)
  }
  ########################################
  digits = max(3, getOption("digits") - 3)
  if(is.null(object2)){
    # stop("The 'anova' function for the sommer package only works to compare mixed models by likelihood ratio tests (LRT), was not intended to provide regular sum of squares output.")
    result <- sequential.fit(object,type=type)
  }else{
    #if(object$maxim){ # user used REML=TRUE, not possible to do LRT
    #  stop("Please fit the models using ML instead of REML by setting the argument REML=FALSE and try again")
    #}else{ #if user used REML=FALSE, then proceed
    if(object$method != object2$method){
      stop("Error! When comparing models please use the same method for the fitted models.")
    }else{
      yu <- summary(object)
      yu2 <- summary(object2)
      dis=c(dim(yu$varcomp)[1]+dim(object$Beta)[1]*dim(object$Beta)[2],
            dim(yu2$varcomp)[1]+dim(object2$Beta)[1]*dim(object2$Beta)[2]) # dimensions
      mods=c("mod1","mod2")
      lls=c( yu$logo[1,1],  yu2$logo[1,1]) # likelihoods
      aics=c(object$AIC, object2$AIC) # AIC's
      bics=c(object$BIC, object2$BIC) # AIC's
      vv=which(dis == max(dis))[1] # which has more variance components BIGGER
      vv2=c(1:2)[which(c(1:2)!= vv)] # SMALLER
      LR = (lls[vv] - lls[vv2])
      r.stat= abs(-2*((LR))) # -2(LL1 - LL2)
      df=dis[vv]-dis[vv2]
      chichi=pchisq((r.stat), df, lower.tail=FALSE)
      if(chichi > 1e-5){
        chichi <- round(chichi,5)
      }
      chichi2=paste(as.character(chichi),signifo(chichi), sep=" ")
      ### construct the table
      cat("Likelihood ratio test for mixed models\n")
      cat("==============================================================\n")
      result=data.frame(Df=c(dis[vv],dis[vv2]), AIC=c(aics[vv],aics[vv2]), 
                        BIC=c(bics[vv],bics[vv2]), loLik=c(lls[vv],lls[vv2]), 
                        Chisq=c("",as.character(round(r.stat,5))), 
                        ChiDf=c("",as.character(df)), PrChisq=c("",chichi2 ))
      rownames(result) <- c(mods[vv],mods[vv2])
      print(result)
      # print(result)
      cat("==============================================================\n")
      cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    }
    #}
  }
  return(result)
}
#### =========== ####
## PLOTING FUNCTION #
#### =========== ####
plot.mmer <- function(x, stnd=TRUE, ...) {
  digits = max(3, getOption("digits") - 3)
  transp <- function (col, alpha = 0.5){
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255,c[3]/255, alpha))
    return(res)
  }
  # provisional model
  prov <- mmer(fixed=x$call$fixed,
               #random=x$call$random,
               rcov=x$call$rcov,
               data=x$data, return.param = TRUE,#reshape.results=TRUE,
               na.method.Y = x$call$na.method.Y, 
               na.method.X = x$call$na.method.X)
  Xlist <- list()
  for(o in 1:length(prov[[3]])){
    Xlist[[o]] <- kronecker(prov[[2]][[o]],prov[[3]][[o]])
  }
  Xm <- do.call(cbind,Xlist)
  # std vs residuals, QQplot (std vs teor quantiles), sqrt(std residuals) vs fitted, std res vs leverage = cook's distance
  traits <- ncol(x$fitted)
  layout(matrix(1:4,2,2))
  for(i in 1:traits){
    plot(x$fitted[,i],scale(x$residuals[,i]),pch=20, col=transp("cadetblue"), ylab="Std Residuals", xlab="Fitted values", main="Residual vs Fitted", bty="n", ...); grid()
    plot(x$fitted[,i],sqrt(abs(scale(x$residuals[,i]))),pch=20, col=transp("thistle4"), ylab="Sqrt Abs Std Residuals", xlab="Fitted values", main="Scale-Location", bty="n", ...); grid()
    qqnorm(scale(x$residuals), pch=20, col=transp("tomato1"), ylab="Std Residuals", bty="n",...); grid()
    hat <- Xm%*%solve(t(Xm)%*%x$Vi%*%Xm)%*%t(Xm)%*%x$Vi # leverage including variance from random effects H= X(X'V-X)X'V-
    plot(diag(hat), scale(x$residuals), pch=20, col=transp("springgreen3"), ylab="Std Residuals", xlab="Leverage", main="Residual vs Leverage", bty="n", ...); grid()
  }
  #####################
  layout(matrix(1,1,1))
}

##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "2.1"))
    stop("This package requires R 2.1 or later")
  assign(".sommer.home", file.path(library, pkg),
         pos=match("package:sommer", search()))
  sommer.version = "4.0 (2019-07-01)" # usually 2 months before it expires
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### check which version is more recent
  #yyy <- 1.8
  #chooseCRANmirror(ind=114)
  #xxx <- available.packages(contriburl = contrib.url(repos="http://mirror.las.iastate.edu/CRAN/", type = getOption("pkgType")))
  
  #xxx <- available.packages()
  #current <- as.numeric(xxx["sommer","Version"])
  ### final check
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  assign(".sommer.version", sommer.version, pos=match("package:sommer", search()))
  if(interactive())
  {
    packageStartupMessage(cyan(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(cyan(paste("[]   Solving Mixed Model Equations in R (sommer) ", sommer.version, "   []",sep="")),appendLF=TRUE)
    packageStartupMessage(cyan(paste("[]   ------------ Multivariate Linear Mixed Models --------------   []")),appendLF=TRUE)
    packageStartupMessage(cyan("[]   Author: Giovanny Covarrubias-Pazaran                           []"),appendLF=TRUE)
    packageStartupMessage(cyan("[]   Published: PLoS ONE 2016, 11(6):1-15                           []"),appendLF=TRUE)
    packageStartupMessage(cyan("[]   Dedicated to the University of Chapingo and the UW-Madison     []"),appendLF=TRUE)
    packageStartupMessage(cyan("[]   Type 'vignette('sommer.start')' for a short tutorial           []"),appendLF=TRUE)
    packageStartupMessage(cyan("[]   Type 'citation('sommer')' to know how to cite sommer           []"),appendLF=TRUE)
    packageStartupMessage(cyan(paste("[]==================================================================[]")),appendLF=TRUE)
    packageStartupMessage(cyan("sommer is updated on CRAN every 4-months due to CRAN policies"),appendLF=TRUE)
    packageStartupMessage(cyan("Newest source is available at https://github.com/covaruber/sommer"),appendLF=TRUE)
    packageStartupMessage(cyan("To install type: library(devtools); install_github('covaruber/sommer')"),appendLF=TRUE)
    
    #if(yyy > current){ # yyy < current in CRAN
    #  packageStartupMessage(paste("Version",current,"is now available."),appendLF=TRUE) # version current
    #  packageStartupMessage(paste("Please update 'sommer' installing the new version."),appendLF=TRUE) # version current
    #}
    #print(image(diag(10),main="sommer 4.0"))
  }
  invisible()
}

#.onLoad <- function(library, pkg) {
#  data(x)
#  library(audio)
#  packageStartupMessage(play(x))
#}

# plot.variogram.mmer <- function(x, stnd=TRUE, ...) {
#   
#   #for(u in 1:length(x)){
#   
#   x0 <- x#[[u]]
#   min.length = 30
#   if(stnd){x0$data$value <- scale(x0$data$value)}
#   values <- matrix(replace(x0$data$value, x0$data$length < min.length,
#                            NA), ncol = length(x0$col.displacement), nrow = length(x0$row.displacement),
#                    byrow = TRUE)
#   
#   print(wireframe(values,drape=TRUE, #main=names(x)[u],
#                   aspect = c(61/87, 0.4), #colorkey=TRUE,
#                   light.source = c(10,0,10), #shade=TRUE,
#                   #screen = list(x0 = -60, y = 50, z=20),
#                   col.regions = topo.colors(100))
#   )
#   # light.source = c(0,0,10),
#   # #region = TRUE,
#   # col.regions = terrain.colors(100),
#   # screen = list(x0 = -60, y = 50, z=20)))
#   
#   # persp(x0, y, z, theta = 135, phi = 30, col = colorRampPalette(c("blue", "pink"))(9500), scale = FALSE,
#   #       ltheta = -120, shade = 0.75, border = NA, box = FALSE)
#   
#   #}
#   # plot3Drgl::persp3Drgl(x0$row.displacement, x0$col.displacement,
#   #                       values, xlab = "Row displacement", ylab = "Col displacement",
#   #                       zlab = "", ticktype = "detailed", col = jet.colors(100))
# }
# 
# variogram.mmer <- function (x, xcoor="ROW", ycoor="RANGE", zcoor=NULL, by=NULL, ...){
#   #x is a dataset with columns Residuals, xcoor, ycoor
#   # xcoor is the name of column for the xcoordinate
#   # ycoor is the name of column for the ycoordinate
#   #library(data.table)
#   
#   x0 <- x$data
#   
#   if(is.null(by)){
#     x0$FIELDINST <- "F1"
#     x0$FIELDINST <- as.factor(x0$FIELDINST)
#     by="FIELDINST"
#   }else{
#     v <- which(colnames(x0) == by)
#     if(length(v)==0){stop("by argument not found in the column names of x0", call. = FALSE)}
#     x0[,by] <- as.factor(x0[,by])
#   }
#   
#   if(is.null(zcoor)){
#     zcoor <- "Residuals"
#     x0[,zcoor] <- x$res.ordim
#   }else{
#     www <- which(names(x$Zus) %in% zcoor)
#     if(length(www)==1){
#       x0[,zcoor] <- x$Zus[[zcoor]]
#     }else if(length(www)>1){
#       zcoor <- paste(zcoor, collapse = "_")
#       gggg <- Reduce("+",x$Zus[www])
#       x0[,zcoor] <- gggg
#     }else{
#       stop("zcoor not found in your model", call. = FALSE)
#     }
#   }
#   
#   are <- which(colnames(x0) %in% c(xcoor,ycoor,zcoor))
#   if(length(are) < 3){
#     stop("One or more of xcoor, ycoor, zcoor don't match your data frame.", call. = FALSE)
#   }
#   
#   ## now add the zcoor
#   
#   x1<- split(x0, f=x0[,by])
#   
#   multires <- lapply(x1, function(x){
#     x.coord <- x[, xcoor]
#     y.coord <- x[, ycoor]
#     residuals <- x[, zcoor]
#     columns <- seq(min(x.coord), max(x.coord), by = min(diff(sort(unique(x.coord)))))
#     rows <- seq(min(y.coord), max(y.coord), by = min(diff(sort(unique(y.coord)))))
#     xy.coord <- data.table(expand.grid(columns = columns, rows = rows))
#     setkeyv(xy.coord, c("rows", "columns"))
#     df <- data.table(columns = x.coord, rows = y.coord, residuals = residuals)
#     setkeyv(df, c("rows", "columns"))
#     df <- df[xy.coord]
#     df <- df[order(df$columns, df$rows), ]
#     resdiff <- c(outer(df$residuals, df$residuals, function(x, 
#                                                             y) 0.5 * (x - y)^2))
#     coldiff <- c(outer(df$columns, df$columns, function(x, y) abs(x - 
#                                                                     y)))
#     coldiff.u <- unique(coldiff)
#     rowdiff <- c(outer(df$rows, df$rows, function(x, y) abs(x - 
#                                                               y)))
#     rowdiff.u <- unique(rowdiff)
#     subsets <- split(resdiff, f = list(coldiff, rowdiff))
#     value <- sapply(subsets, mean, na.rm = TRUE)
#     length <- sapply(subsets, function(x) sum(!is.na(x)))
#     length[-1] <- length[-1]/2
#     res <- list(data = data.frame(value = value, length = length), 
#                 col.displacement = coldiff.u, row.displacement = rowdiff.u)
#     class(res) <- "variogram.mmer"
#     return(res)
#   })
#   #print(str(multires))
#   # multires <- as.data.frame(do.call(rbind,multires))
#   # class(multires) <- "variogram.sommer"
#   # 
#   # funny <- as.formula(paste(zcoor,"~",xcoor,"*",ycoor,"|",by, sep="" ))
#   # multires2 <- list(data=multires,funny=funny)
#   class(multires) <- "variogram.mmer"
#   
#   return(multires)
# }

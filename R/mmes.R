mmes <- function(fixed, random, rcov, data, W,
                 nIters=50, tolParConvLL = 1e-04,
                 tolParConvNorm = 1e-04, tolParInv = 1e-06,
                 naMethodX="exclude",
                 naMethodY="exclude",
                 returnParam=FALSE,
                 dateWarning=TRUE,
                 verbose=TRUE, 
                 addScaleParam=NULL,
                 stepWeight=NULL, emWeight=NULL, 
                 contrasts=NULL,
                 getPEV=TRUE, henderson=FALSE){
  
  my.date <- "2025-08-01"
  your.date <- Sys.Date()
  ## if your month is greater than my month you are outdated
  if(dateWarning){
    if (your.date > my.date) {
      cat("Version out of date. Please update sommer to the newest version using:\ninstall.packages('sommer') in a new session\n Use the 'dateWarning' argument to disable the warning message.")
    }
  }
  
  if(missing(data)){
    data <- environment(fixed)
    if(!missing(random)){
      data2 <- environment(random)
    }
    nodata <-TRUE
    cat("data argument not provided \n")
  }else{nodata=FALSE; data <- as.data.frame(data)}
  
  if(missing(rcov)){
    rcov = as.formula("~units")
  }
  
  #################
  ## do the needed for naMethodY and naMethodX
  dataor <- data
  provdat <- subdata(data, fixed=fixed, na.method.Y = naMethodY,na.method.X=naMethodX)
  data <- provdat$datar
  nonMissing <- provdat$good
  #################
  data$units <- levels(as.factor(paste("u",1:nrow(data),sep="")))
  #################
  ## get Y matrix
  response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
  responsef <- as.formula(paste(response,"~1"))
  mfna <- try(model.frame(responsef, data = data, na.action = na.pass), silent = TRUE)
  if (is(mfna, "try-error") ) { # class(mfna) == "try-error"
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", call. = FALSE)
  }
  mfna <- eval(mfna, data, parent.frame())
  yvar <- sparse.model.matrix(as.formula(paste("~",response,"-1")),data)
  nt <- ncol(yvar)
  if(nt==1){colnames(yvar) <- response}
  Vy <- var(yvar[,1])
  # yvar <- scale(yvar)
  #################
  ## get Zs and Ks
  
  Z <- Ai <- theta <- thetaC <- thetaF <- sp <- list()
  Zind <- numeric()
  rTermsNames <- list()
  counter <- 1
  if(!missing(random)){ # if there's random effects
    
    yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]] # random parts
    rtermss <- apply(data.frame(yuyu),1,function(x){ # split random terms
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    # print(rtermss)
    for(u in 1:length(rtermss)){ # for each random effect u=1
      checkvs <- intersect(all.names(as.formula(paste0("~",rtermss[u]))),c("vsm","spl2Dc")) # which(all.names(as.formula(paste0("~",rtermss[u]))) %in% c("vs","spl2Da","spl2Db")) # grep("vs\\(",rtermss[u])
      
      if(length(checkvs)==0){ ## if this term is not in a variance structure put it inside
        rtermss[u] <- paste("sommer::vsm( sommer::ism(",rtermss[u],") )")
      }
      ff <- eval(parse(text = rtermss[u]),data,parent.frame()) # evaluate the variance structure
      Z <- c(Z, lapply(ff$Z, function(x){if(nrow(x) != length(nonMissing)){return(x[nonMissing,])}else{return(x)} }) )
      # Z <- c(Z, ff$Z)
      Ai <- c(Ai, ff$Gu)
      theta[[u]] <- ff$theta
      thetaC[[u]] <- ff$thetaC
      thetaF[[u]] <- ff$thetaF
      sp[[u]] <- ff$sp # rep(ff$sp,length(which(ff$thetaC > 0)))
      Zind <- c(Zind, rep(u,length(ff$Z)))
      checkvs <- numeric() # restart the check
      ## names for monitor
      baseNames <- which( ff$thetaC > 0, arr.ind = TRUE)
      s1 <- paste(rownames(ff$thetaC)[baseNames[,"row"]], colnames(ff$thetaC)[baseNames[,"col"]],sep = ":")
      s2 <- paste(all.vars(as.formula(paste("~",rtermss[u]))),collapse=":")
      rTermsNames[[u]] <- paste(s2,s1,sep=":")
      
      counter <- counter + 1
    }
  }
  
  # nEffects <- sum(unlist(lapply(Z,ncol)))
  # nRecords <- length(yvar)
  #################
  ## get Rs
  
  yuyur <- strsplit(as.character(rcov[2]), split = "[+]")[[1]]
  rcovtermss <- apply(data.frame(yuyur),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  
  S <- list()
  Spartitions <- list()
  for(u in 1:length(rcovtermss)){ # for each random effect
    checkvs <- intersect(all.names(as.formula(paste0("~",rcovtermss[u]))),c("vsm","gvs","spl2Da","spl2Db")) # which(all.names(as.formula(paste0("~",rtermss[u]))) %in% c("vs","spl2Da","spl2Db")) # grep("vs\\(",rtermss[u])
    
    if(length(checkvs)==0){ ## if this term is not in a variance structure put it inside
      rcovtermss[u] <- paste("sommer::vsm( sommer::ism(",rcovtermss[u],") )")
    }
    
    ff <- eval(parse(text = rcovtermss[u]),data,parent.frame()) # evalaute the variance structure
    S <- c(S, ff$Z)
    Spartitions <- c(Spartitions, ff$partitionsR)
    ## constraint
    residualsNonFixed <- which(ff$thetaC != 3, arr.ind = TRUE)
    if(nrow(residualsNonFixed) > 0){
      ff$theta[residualsNonFixed] <- ff$theta[residualsNonFixed] * 5
    }
    theta[[counter]] <- ff$theta
    thetaC[[counter]] <- ff$thetaC
    thetaF[[counter]] <- ff$thetaF
    sp[[counter]] <- ff$sp#rep(ff$sp,length(ff$Z))
    
    baseNames <- which( ff$thetaC > 0, arr.ind = TRUE)
    s1 <- paste(rownames(ff$thetaC)[baseNames[,"row"]], colnames(ff$thetaC)[baseNames[,"col"]],sep = ":")
    s2 <- paste(all.vars(as.formula(paste("~",rcovtermss[u]))),collapse=":")
    rTermsNames[[counter]] <- paste(s2,s1,sep=":")
    
    checkvs <- numeric() # restart the check
    counter <- counter + 1
  }
  #################
  #################
  ## get Xs
  data$`1` <-1
  newfixed=fixed
  fixedTerms <- gsub(" ", "", strsplit(as.character(fixed[3]), split = "[+-]")[[1]])
  mf <- try(model.frame(newfixed, data = data, na.action = na.pass), silent = TRUE)
  mf <- eval(mf, parent.frame())
  X <-  Matrix::sparse.model.matrix(newfixed, mf, contrasts.arg=contrasts)
  
  
  partitionsX <- list()#as.data.frame(matrix(NA,length(fixedTerms),2))
  for(ix in 1:length(fixedTerms)){ # save indices for partitions of each fixed effect
    effs <- colnames(Matrix::sparse.model.matrix(as.formula(paste("~",fixedTerms[ix],"-1")), mf))
    effs2 <- colnames(Matrix::sparse.model.matrix(as.formula(paste("~",fixedTerms[ix])), mf))
    partitionsX[[ix]] <- matrix(which(colnames(X) %in% c(effs,effs2)),nrow=1)
  }
  names(partitionsX) <- fixedTerms
  classColumns <- lapply(data,class)
  
  for(ix in 1:length(fixedTerms)){ # clean column names in X matrix
    colnamesBase <- colnames(X)[partitionsX[[ix]]]
    colnamesBaseList <- strsplit(colnamesBase,":")
    toRemoveList <- strsplit(fixedTerms[ix],":")[[1]] # words to remove from the level names in the ix.th fixed effect
    # print(toRemoveList)
    if("1" %in% unlist(toRemoveList)){}else{toRemoveList <- all.vars(as.formula(paste("~",paste(toRemoveList, collapse = "+"))))}
    for(j in 1:length(toRemoveList)){
      if( toRemoveList[[j]] %in% names(classColumns) ){
        
        if( classColumns[[toRemoveList[[j]]]] != "numeric"){ # only remove the name from the level if is structure between factors, not for random regressions
          nc <- nchar(gsub(" ", "", toRemoveList[[j]], fixed = TRUE)) # number of letters to remove
          colnamesBaseList <- lapply(colnamesBaseList, function(h){
            if(is.na(h[j])){ # is the intercept? no
              return(h)
            }else{ # is the intercept? yes
              if(length(grep(toRemoveList[[j]],h[j])) == 1){ # if the factor word matches in the level
                if(nchar(h[j]) > nc){h[j] <- substr(h[j],1+nc,nchar(h[j]))}
              }
              return(h)
            }
          }) # only remove the initial name if the name is actually longer
        }
        
      }
    }
    colnames(X)[partitionsX[[ix]]] <- unlist(lapply(lapply(colnamesBaseList,na.omit), function(x){paste(x, collapse=":")}))
  }
  step1 <- gsub(" ", "", strsplit(as.character(fixed[3]), split = "[-]")[[1]])
  step2 <- unlist(apply(data.frame(step1),1,function(x){strsplit(as.character(x), split = "[+]")[[1]]}))
  intercCheck <- ifelse(length(intersect(c("1","-1"),step2)) == 0, TRUE, FALSE) # if length is zero it means that we have an intercept
  if(intercCheck){colnames(X)[1] <- "Intercept"}
  
  #################
  #################
  ## weight matrix
  
  if(missing(W)){ # provide W but don't use it
    x <- data.frame(d=as.factor(1:length(yvar)))
    W <- sparse.model.matrix(~d-1, x)
    useH=FALSE
  }else{
    W <- as(as(as( W ,  "dMatrix"), "generalMatrix"), "CsparseMatrix") # as(W, Class = "dgCMatrix")
    useH=TRUE
  }
  
  #################
  #################
  ## information weights
  
  if(is.null(emWeight)){
    if(henderson==FALSE){ # p > n
      emWeight <- rep(0, nIters)
    }else{ # n > p
      # initialEmSteps <- logspace(round(nIters*.8), 1, 0.009) # 80% of the iterations requested are used for the logarithmic decrease
      # restEmSteps <- rep(0.009, nIters - length(initialEmSteps)) # the rest we assign a very small emWeight value
      emWeight <- stan(logspace(seq(1,-1,- 2/nIters), p=3)) #c( initialEmSteps, restEmSteps) # plot(emWeight) # we bind both for the modeling
    }
    
  }
  if(is.null(stepWeight)){
    w <- which(emWeight <= .5) # where AI starts
    if(length(w) > 1){ # w has at least length = 2
      stepWeight <- rep(.9,nIters);
      if(nIters > 1){stepWeight[w[1:2]] <- c(0.5,0.7)} # .5, .7
    }else{
      stepWeight <- rep(.9,nIters);
      if(nIters > 1){stepWeight[1:2] <- c(0.5,0.7)} # .5, .7
    }
  }
  
  #################
  #################
  ## information weights
  theta <- lapply(theta, function(x){return(x*Vy)})
  
  thetaFinput <- do.call(adiag1,thetaF)
  if(is.null(addScaleParam)){addScaleParam=0}
  thetaFinputSP <- unlist(sp)
  thetaFinput <- cbind(thetaFinput,thetaFinputSP)
  thetaFinput
  
  if(henderson){
    nInverses <- length(which(unlist(lapply(Ai, function(x){attributes(x)$inverse})) == TRUE))
    if(nInverses != length(Ai)){
      stop("You have selected the 'henderson' algorithm which requires all relationship
      matrices to be inverted. Please make sure that you have inverted your
      matrices and set the attributes of your matrices as follows:
           Gu = as(as(as( Gu,  'dMatrix'), 'generalMatrix'), 'CsparseMatrix')
           attr(Gu, 'inverse')=TRUE 
      where 'Gu' is to be replaced with the name of your matrix.", call. = FALSE)
    }
  }
  # else{
  #   nNoInverses <- length(which(unlist(lapply(Ai, function(x){attributes(x)$inverse})) == FALSE))
  #   if(nNoInverses != length(Ai)){
  #     stop("You have selected the 'direct-inversion' algorithm which requires all 
  #     relationship matrices to NOT be inverted. Please make sure that you have 
  #     provided your raw matrices and set the attributes of your matrices as follows:
  #          Gu = as(as(as( Gu,  'dMatrix'), 'generalMatrix'), 'CsparseMatrix')
  #          attr(Gu, 'inverse')=FALSE 
  #     where 'Gu' is to be replaced with the name of your matrix.", call. = FALSE)
  #   }
  # }
  
  if(returnParam){ # if user just wants to get input matrices
    
    res <- list(yvar=yvar, X=X,Z=Z,Zind=Zind,Ai=Ai,S=S,Spartitions=Spartitions, W=W, useH=useH,
                nIters=nIters, tolParConvLL=tolParConvLL, tolParConvNorm=tolParConvNorm,
                tolParInv=tolParInv,
                verbose=verbose, addScaleParam=addScaleParam,
                theta=theta,thetaC=thetaC, thetaFinput=thetaFinput,
                stepWeight=stepWeight,emWeight=emWeight, 
                rtermss=rtermss, partitionsX=partitionsX, getPEV=getPEV, rTermsNames=rTermsNames
    )
    
    # yvar=res$yvar; X=res$X;Z=res$Z;Zind=res$Zind;Ai=res$Ai;S=res$S;
    # Spartitions=res$Spartitions; W=res$W; useH=res$useH;
    # nIters=res$nIters; tolParConvLL=res$tolParConvLL; tolParConvNorm=res$tolParConvNorm;
    # tolParInv=res$tolParInv;
    # verbose=res$verbose; addScaleParam=res$addScaleParam;
    # theta=res$theta;thetaC=res$thetaC; thetaFinput=res$thetaFinput;
    # stepWeight=res$stepWeight;emWeight=res$emWeight; 
    # rtermss=res$rtermss; partitionsX=res$partitionsX; getPEV=res$getPEV
    
  }else{
    
    #################################################
    if(henderson == FALSE){ # p > n = DIRECT INVERSION: transform all matrices to expected format
      # get isInvW
      isInvW=FALSE
      AI=FALSE # use newton raphson
      returnScaled=FALSE # return scaled variance parameters
      # translate vsm S into vsm R
      R <- rep(list(Matrix::Diagonal(x= rep(0, nrow(yvar)) )), length(S) )
      for(iR in 1:length(S)){ # iR=1
        R[[iR]][Spartitions[[iR]][1,1]:Spartitions[[iR]][1,2],
                Spartitions[[iR]][1,1]:Spartitions[[iR]][1,2] ] = S[[iR]]
      }
      R <- lapply(R,function(x){as(as(as( x,  "dMatrix"), "generalMatrix"), "CsparseMatrix")})
      # translate vsm Z into vsm Z
      
      THETA <- THETAc <- K <- Zdi <- list(); counter=1
      vary <- var(yvar[,1])
      thetaIndex <- thetaConstIndex <- numeric()
      
      
      for(iTheta in 1:length(theta)){
        for(iRow in 1:nrow(theta[[iTheta]])){ # iRow=1
          for(iCol in iRow:ncol(theta[[iTheta]])){ # iCol=1
            if(theta[[iTheta]][iRow,iCol] != 0){
              THETA[[counter]] <- matrix(theta[[iTheta]][iRow,iCol]/vary,1,1)
              THETAc[[counter]] <- matrix(thetaC[[iTheta]][iRow,iCol],1,1)
              thetaConstIndex[counter] <- thetaC[[iTheta]][iRow,iCol]
              thetaIndex[counter] <- iTheta
              colnames(THETA[[counter]] ) <- rownames(THETA[[counter]] ) <- colnames(THETAc[[counter]] ) <- rownames(THETAc[[counter]] ) <- paste( colnames(theta[[iTheta]])[iRow], colnames(theta[[iTheta]])[iCol], sep="_" )
              # if variance components only (we don't do this for residual vcs)
              if(iTheta < length(theta)){ 
                useZs <- which(Zind == iTheta)
                if(iRow == iCol){ # if variance component
                  K[[counter]] <- Ai[[iTheta]]
                  Zdi[[counter]] <-  Z[[useZs[iRow]]] 
                }else{ # if covariance component
                  nrowcol <- nrow(Ai[[iTheta]])
                  Kcov <- Matrix::Matrix(0,nrowcol*2,nrowcol*2)
                  Kcov <- as(as(as( Kcov,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
                  Kcov[1:nrowcol,(nrowcol+1):nrow(Kcov)] = Ai[[iTheta]]
                  Kcov[(nrowcol+1):nrow(Kcov),1:nrowcol] = Ai[[iTheta]]
                  K[[counter]] <- Kcov
                  Zdi[[counter]] <-  cbind( Z[[useZs[iRow]]], Z[[useZs[iCol]]] ) 
                }
              }
              # enf of if variance component
              counter=counter+1
            }
          }
        }
      }
      
      if(!missing(random)){
        XZ <- cbind(X,do.call(cbind,Z))
      }else{XZ <- X}
      theta <- THETA; THETA <- NULL
      
      res <- .Call("_sommer_newton_di_sp",PACKAGE = "sommer",
                   yvar,
                   list(X),
                   list(matrix(1)),
                   Zdi,K,R,
                   theta,THETAc,
                   W,
                   isInvW,
                   nIters, tolParConvLL, tolParInv,
                   AI,getPEV,verbose, returnScaled,
                   stepWeight, emWeight,
                   thetaC, thetaIndex)
      
      # res <- newton_di_sp(
      #              yvar,
      #              list(X),
      #              list(matrix(1)),
      #              Zdi,K,R,
      #              theta,THETAc,
      #              W,
      #              isInvW,
      #              nIters, tolParConvLL, tolParInv,
      #              AI,getPEV,verbose, returnScaled,
      #              stepWeight, emWeight,
      #              thetaC, thetaIndex)
      
    }else if(henderson == TRUE){ # n > p HENDERSON
      
      Si <- lapply(S, solve)
      res <- .Call("_sommer_ai_mme_sp",PACKAGE = "sommer",
                   X,Z, Zind,
                   Ai,yvar,
                   S, Spartitions, W, useH,
                   nIters, tolParConvLL, tolParConvNorm,
                   tolParInv,theta,
                   thetaC,thetaFinput,
                   addScaleParam,
                   emWeight,
                   stepWeight,
                   verbose)
      
      # res <- ai_mme_sp(
      #              X,Z, Zind,
      #              Ai,yvar,
      #              Si, Spartitions, W, useH,
      #              nIters, tolParConvLL, tolParConvNorm,
      #              tolParInv,theta,
      #              thetaC,thetaFinput,
      #              addScaleParam,
      #              emWeight,
      #              stepWeight,
      #              verbose)
      
    }
    ###### add rownames and build uList
    rownames(res$b) <- colnames(X)
    res$thetaC <- thetaC # also residuals
    if(!missing(random)){
      rownames(res$u) <- unlist(lapply(Z, colnames))
    }
    rownames(res$bu) <- c(rownames(res$b),rownames(res$u))
    rownames(res$monitor) <- unlist(rTermsNames)
    res$data <- data
    res$y <- yvar
    res$partitionsX <- partitionsX

    if(!missing(random)){ # mock
      names(res$theta) <- names(res$thetaC) <- c(rtermss,"units")
      
      names(res$partitions) <- rtermss
      names(res$uList) <- names(res$uPevList) <- rtermss
      for(i in 1:length(res$partitions)){ # i=1
        colnames(res$uList[[i]]) <- colnames(thetaC[[i]])
        rownames(res$uList[[i]]) <- rownames(res$bu)[res$partitions[[i]][1,1]:res$partitions[[i]][1,2]]
        if(getPEV){
          colnames(res$uPevList[[i]]) <- colnames(thetaC[[i]])
          rownames(res$uPevList[[i]]) <- rownames(res$bu)[res$partitions[[i]][1,1]:res$partitions[[i]][1,2]]
        }
      }
      
      if(henderson==FALSE){ ######## adding ulist and upevlist similar to henderson mmes
        res$W <- XZ
      }

    } # enf of 'if(missing(random))'
    ## adding D table for predictions
    if(!missing(random)){
      res$args <- list(fixed=fixed, random=random, rcov=rcov)
      res$Dtable <- data.frame(type=c(rep("fixed",length(res$partitionsX)),
                                      rep("random",length(res$partitions))
                                      ),
                               term=c(names(res$partitionsX),names(res$partitions)),
                               include=FALSE,average=FALSE)
    }else{
      res$args <- list(fixed=fixed, rcov=rcov)
      res$Dtable <- data.frame(type=c(rep("fixed",length(res$partitionsX))),term=c(names(res$partitionsX),names(res$partitions)),include=FALSE,average=FALSE)
    }
    # 
    class(res)<-c("mmes")
  }
  return(res)
}

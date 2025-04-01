mme <- function(fixed, random, rcov, data, W,
                 nIters=25, tolParConvLL = 1e-03,
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
  
  my.date <- "2025-06-01"
  your.date <- Sys.Date()
  ## if your month is greater than my month you are outdated
  if(dateWarning){
    if (your.date > my.date) {
      cat("Version out of date. Please update sommex to the newest version using:\ninstall.packages('sommex') in a new session\n Use the 'dateWarning' argument to disable the warning message.")
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
        rtermss[u] <- paste("vsm( ism(",rtermss[u],") )")
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
      rcovtermss[u] <- paste("vsm( ism(",rcovtermss[u],") )")
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
    if(mme==FALSE){ # p > n
      emWeight <- rep(0, nIters)
    }else{ # n > p
      initialEmSteps <- logspace(round(nIters*.8), 1, 0.009) # 80% of the iterations requested are used for the logarithmic decrease
      restEmSteps <- rep(0.009, nIters - length(initialEmSteps)) # the rest we assign a very small emWeight value
      emWeight <- c( initialEmSteps, restEmSteps) # plot(emWeight) # we bind both for the modeling
    }
    
  }
  if(is.null(stepWeight)){
    w <- which(emWeight <= .04) # where AI starts
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
  
  if(mme){
    nInverses <- length(unlist(lapply(Ai, function(x){attributes(x)$inverse})))
    if(nInverses != length(Ai)){
      stop("You have selected the 'mme' algorithm which requires all relationship
      matrices to be inverted. Please make sure that you have inverted your
      matrices and set the attribute to your matrices as follows:
           attr(Gu, 'inverse')=TRUE 
      where 'Gu' is to be replaced with your matrix name.", call. = FALSE)
    }
  }
  
  if(returnParam){ # if user just wants to get input matrices
    
    res <- list(yvar=yvar, X=X,Z=Z,Zind=Zind,Ai=Ai,S=S,Spartitions=Spartitions, W=W, useH=useH,
                nIters=nIters, tolParConvLL=tolParConvLL, tolParConvNorm=tolParConvNorm,
                tolParInv=tolParInv,
                verbose=verbose, addScaleParam=addScaleParam,
                theta=theta,thetaC=thetaC, thetaF=thetaFinput,
                stepWeight=stepWeight,emWeight=emWeight
    )
    
  }else{
    
    #################################################
    if(mme == FALSE){ # p > n = DIRECT INVERSION: transform all matrices to expected format
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
                  K[[counter]] <- as.matrix(Ai[[iTheta]])
                  Zdi[[counter]] <-  Z[[useZs[iRow]]] 
                }else{ # if covariance component
                  nrowcol <- nrow(Ai[[iTheta]])
                  Kcov <- matrix(0,nrowcol*2,nrowcol*2)
                  Kcov[1:nrowcol,(nrowcol+1):nrow(Kcov)] = as.matrix(Ai[[iTheta]])
                  Kcov[(nrowcol+1):nrow(Kcov),1:nrowcol] = as.matrix(Ai[[iTheta]])
                  K[[counter]] <- as.matrix(Kcov)
                  Zdi[[counter]] <-  cbind( Z[[useZs[iRow]]], Z[[useZs[iCol]]] ) 
                }
              }
              # enf of if variance component
              counter=counter+1
            }
          }
        }
      }
      # print(thetaIndex)
      
      if(!missing(random)){
        # print(length(nonMissing))
        XZY <- cbind(X,do.call(cbind,Z))
      }else{XZY <- NULL}
      Z <- Zdi; Zdi <- NULL
      theta <- THETA; THETA <- NULL
      # thetaC <- THETAc;  THETAc <- NULL
      
      res <- .Call("_sommex_MNR",PACKAGE = "sommex",
                   as.matrix(yvar), 
                   list(as.matrix(X)),
                   list(matrix(1)),
                   Z,K,R,
                   theta,THETAc, 
                   as.matrix(W), 
                   isInvW,
                   nIters, tolParConvLL, tolParInv,
                   AI,getPEV,verbose, returnScaled, stepWeight, emWeight)
      
    }else if(mme == TRUE){ # n > p HENDERSON
      
      res <- .Call("_sommex_ai_mme_sp",PACKAGE = "sommex",
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
      
    }
    ###### add rownames and build uList
    rownames(res$b) <- colnames(X)
    if(!missing(random)){
      rownames(res$u) <- unlist(lapply(Z, colnames))
      res$thetaC <- thetaC
    }
    rownames(res$bu) <- c(rownames(res$b),rownames(res$u))
    # print(unlist(rTermsNames))
    rownames(res$monitor) <- unlist(rTermsNames)
    res$sigma <- res$monitor[,which(res$llik[1,] == max(res$llik[1,]))] # we return the ones with max llik
    res$data <- data
    res$y <- yvar
    res$partitionsX <- partitionsX
    uList <- uPevList <- vector(mode="list",length = length(thetaC)-1)

    if(!missing(random)){
      names(uList) <- names(uPevList) <- rtermss
      if(mme==FALSE){ ######## adding ulist and upevlist similar to mme mme
        names(res$partitions) <- rtermss
        newtheta <- list()
        for(iTheta in 1:(length(thetaC)-1)){ # iTheta=2 # each element in the list
          effsToUse <- which(thetaIndex == iTheta)
          nEffs <- min(unlist(lapply(res$uList0[effsToUse], nrow)))
          blupTable <- pevTable <-  Matrix::Matrix(0, nrow=nEffs, ncol=ncol(thetaC[[iTheta]]))
          newthetasub <- thetaC[[iTheta]]
          counter=1
          for(iRow in 1:nrow(thetaC[[iTheta]])){ # iRow=1
            for(iCol in iRow:ncol(thetaC[[iTheta]])){ # iCol=2

              if(thetaC[[iTheta]][iRow,iCol] != 0){ # if is to be estimated
                newthetasub[iRow,iCol] <- newthetasub[iCol,iRow] <- res$theta[[effsToUse[counter]]]
                if(iRow==iCol){ # variance component

                  blupTable[,iRow] <- blupTable[,iRow] + res$uList0[[effsToUse[counter]]]
                  pevTable[,iRow] <- pevTable[,iRow] + diag(res$Ci[[effsToUse[counter]]])
                  counter=counter+1

                }else{ # covariance component

                  prov <- res$uList0[[effsToUse[counter]]]
                  indexprov <- c( rep(iRow, nrow(prov)/2), rep(iCol, nrow(prov)/2) )
                  # iRow
                  blupTable[,iRow] <- blupTable[,iRow] +  res$uList0[[effsToUse[counter]]][which(indexprov == iRow),]
                  pevTable[,iRow] <- pevTable[,iRow] +  diag(res$Ci[[effsToUse[counter]]])[which(indexprov == iRow)]
                  # iCol
                  blupTable[,iCol] <- blupTable[,iCol] +  res$uList0[[effsToUse[counter]]][which(indexprov == iCol),]
                  pevTable[,iCol] <- pevTable[,iCol] +  diag(res$Ci[[effsToUse[counter]]])[which(indexprov == iCol)]
                  counter=counter+1

                }
              }
            }
          }
          colnames(blupTable) <- colnames(pevTable) <- colnames(thetaC[[iTheta]])
          uList[[iTheta]] <- blupTable
          uPevList[[iTheta]] <- pevTable
          newtheta[[iTheta]] <- newthetasub
        }
        # move the PEV to the Cii matrix
        if(length(Z) > 0){ # there's random effects
          res$Ci <- sommex::adiag1(res$Ci_11, do.call(sommex::adiag1, res$Ci))
          res$W <- XZY
        }else{
          res$Ci <- res$Ci_11
        }
        residualthetas <- which(thetaIndex == length(thetaC))
        newtheta[[iTheta+1]] <- res$theta[residualthetas]
        res$theta <- newtheta
        res$uList0 <- NULL
      }else if(mme==TRUE){ ######## adding ulist and upevlist in mme mme
        names(res$partitions) <- rtermss
        for(i in 1:length(res$partitions)){ # i=1
          blupTable <- apply(res$partitions[[i]],1,function(x2){
            toreturn <- res$bu[(x2[1]):(x2[2]),,drop=FALSE]
            return(toreturn)
          })
          rownames(blupTable) <- rownames(res$bu)[res$partitions[[i]][1,1]:res$partitions[[i]][1,2]]
          pevTable <- apply(res$partitions[[i]],1,function(x2){return(diag(res$Ci)[(x2[1]):(x2[2])])})
          colnames(blupTable) <- colnames(pevTable) <- colnames(theta[[i]])
          rownames(pevTable) <- rownames(blupTable)
          uList[[i]] <- blupTable; uPevList[[i]] <- pevTable
        };
        blupTable=NULL; pevTable=NULL;
      }
      names(res$theta) <- names(res$thetaC) <- c(rtermss,"units")

    }
    res$uList <- uList; res$uPevList <- uPevList
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
    class(res)<-c("mme")
  }
  return(res)
}

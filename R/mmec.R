mmec <- function(fixed, random, rcov, data, W,
                  nIters=20, tolParConvLL = 1e-03, 
                  tolParConvNorm = 0.05, tolParInv = 1e-06,
                  naMethodX="exclude",
                  naMethodY="exclude",
                  returnParam=FALSE,
                  dateWarning=TRUE,
                  verbose=TRUE, addScaleParam=NULL,
                  stepWeight=NULL, emWeight=NULL){
  
  my.year <- 2022
  my.month <- 6 #month when the user will start to get notifications the 1st day of next month
  ### if my month = 5, user will start to get notification in June 1st (next month)
  datee <- Sys.Date()
  year.mo.day <- as.numeric(strsplit(as.character(datee),"-")[[1]])# <- as.numeric(strsplit(gsub("....-","",datee),"-")[[1]])
  your.year <- year.mo.day[1]
  your.month <- year.mo.day[2]
  ## if your month is greater than my month you are outdated
  if(dateWarning){
    if(your.month > my.month & your.year >= my.year){
      # error if your month is greater and your year is smaller
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
  }else{nodata=FALSE}
  
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
  
  Z <- Ai <- theta <- thetaC <- thetaF <- list()
  Zind <- numeric()
  rTermsNames <- list()
  counter <- 1
  if(!missing(random)){ # if there's random effects
    
    yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]] # random parts
    rtermss <- apply(data.frame(yuyu),1,function(x){ # split random terms
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    
    for(u in 1:length(rtermss)){ # for each random effect
      checkvs <- intersect(all.names(as.formula(paste0("~",rtermss[u]))),c("vsc","spl2Dc1")) # which(all.names(as.formula(paste0("~",rtermss[u]))) %in% c("vs","spl2Da","spl2Db")) # grep("vs\\(",rtermss[u])
      
      if(length(checkvs)==0){ ## if this term is not in a variance structure put it inside
        rtermss[u] <- paste("vsc( isc(",rtermss[u],") )")
      }
      ff <- eval(parse(text = rtermss[u]),data,parent.frame()) # evaluate the variance structure
      Z <- c(Z, ff$Z)
      Ai <- c(Ai, ff$Gu)
      theta[[u]] <- ff$theta
      thetaC[[u]] <- ff$thetaC
      thetaF[[u]] <- ff$thetaF
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
  
  #################
  ## get Rs
  
  yuyur <- strsplit(as.character(rcov[2]), split = "[+]")[[1]]
  rcovtermss <- apply(data.frame(yuyur),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  
  S <- list()
  Spartitions <- list()
  for(u in 1:length(rcovtermss)){ # for each random effect
    checkvs <- intersect(all.names(as.formula(paste0("~",rcovtermss[u]))),c("vsc","gvs","spl2Da","spl2Db")) # which(all.names(as.formula(paste0("~",rtermss[u]))) %in% c("vs","spl2Da","spl2Db")) # grep("vs\\(",rtermss[u])
    
    if(length(checkvs)==0){ ## if this term is not in a variance structure put it inside
      rcovtermss[u] <- paste("vsc( isc(",rcovtermss[u],") )")
    }
    
    ff <- eval(parse(text = rcovtermss[u]),data,parent.frame()) # evalaute the variance structure
    S <- c(S, ff$Z)
    Spartitions <- c(Spartitions, ff$partitionsR)
    theta[[counter]] <- ff$theta*5
    thetaC[[counter]] <- ff$thetaC
    thetaF[[counter]] <- ff$thetaF
    
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
  
  # Expand the fixed terms to handle interactions
  # fixedtermss <- attr(terms(fixed),"term.labels")
  # if(length(fixedtermss) == 0){fixedtermss="1"}
  # newfixed <- as.formula(paste("~",paste(fixedtermss,collapse = "+")))
  newfixed=fixed
  # fixedTerms <- all.vars(newfixed)[-1]
  fixedTerms <- strsplit(as.character(fixed[3]), split = "[+]")[[1]]
  mf <- try(model.frame(newfixed, data = data, na.action = na.pass), silent = TRUE)
  mf <- eval(mf, parent.frame())
  X <-  Matrix::sparse.model.matrix(newfixed, mf)
  
  partitionsX <- list()#as.data.frame(matrix(NA,length(fixedTerms),2))
  for(ix in 1:length(fixedTerms)){
    effs <- colnames(Matrix::sparse.model.matrix(as.formula(paste("~",fixedTerms[ix],"-1")), mf))
    effs2 <- colnames(Matrix::sparse.model.matrix(as.formula(paste("~",fixedTerms[ix])), mf))
    partitionsX[[ix]] <- matrix(which(colnames(X) %in% c(effs,effs2)),nrow=1)
  }
  names(partitionsX) <- fixedTerms
  #################
  #################
  ## weight matrix
  
  if(missing(W)){
    x <- data.frame(d=as.factor(1:length(yvar)))
    W <- sparse.model.matrix(~d-1, x)
    useH=FALSE
  }else{
    W <- as(W, Class = "sparseMatrix")
    useH=TRUE
  }
  
  #################
  #################
  ## information weights
  
  if(is.null(emWeight)){
    emWeight <- c(seq(1,.1,-.15),rep(0.04,nIters))
  }
  if(is.null(stepWeight)){
    w <- which(emWeight <= .04) # where AI starts
    if(length(w) > 0){
      stepWeight <- rep(.9,nIters); 
      if(nIters > 1){stepWeight[w[1:2]] <- c(0.5,0.7)} # .5, .7
    }else{
      stepWeight <- rep(.9,nIters); 
      if(nIters > 1){stepWeight[1:2] <- c(0.5,0.7)} # .5, .7
    }
    # stepWeight[1:3]=3
  }
  
  #################
  #################
  ## information weights
  theta <- lapply(theta, function(x){return(x*Vy)})

  if(returnParam){
    
    # args <- list(fixed=fixed, random=random, rcov=rcov, data=data, W=W,
    #              nIters=nIters, tolParConv=tolParConv, tolParInv=tolParInv,
    #              naMethodX=naMethodX,
    #              naMethodY=naMethodY,
    #              returnParam=returnParam,
    #              dateWarning=dateWarning,
    #              verbose=verbose)
    # 
    # good <- provdat$good
    
    thetaFinput <- do.call(adiag1,thetaF)
    if(is.null(addScaleParam)){addScaleParam=0}
    thetaFinput <- cbind(thetaFinput,matrix(0,nrow(thetaFinput),length(addScaleParam)))
    thetaFinput
    
    res <- list(yvar=yvar, X=X,Z=Z,Zind=Zind,Ai=Ai,S=S,Spartitions=Spartitions, W=W, useH=useH,
                nIters=nIters, tolParConvLL=tolParConvLL, tolParConvNorm=tolParConvNorm, 
                tolParInv=tolParInv,
                verbose=verbose, addScaleParam=addScaleParam,
                theta=theta,thetaC=thetaC, thetaF=thetaFinput,
                stepWeight=stepWeight,emWeight=emWeight 
    )
  }else{
    
    thetaFinput <- do.call(adiag1,thetaF)
    if(is.null(addScaleParam)){addScaleParam=0}
    thetaFinput <- cbind(thetaFinput,matrix(0,nrow(thetaFinput),length(addScaleParam)))
    thetaFinput
    
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
    
    # res <-ai_mme_sp(X=X,ZI=Z, Zind=Zind,
    #                 AiI=Ai,y=yvar,
    #                 SI=S, Spartitions,
    #                 H=W, useH=useH,
    #                 nIters=nIters, tolParConvLL=tolParConvLL,
    #                 tolParConvNorm=tolParConvNorm,
    #                 tolParInv=tolParInv,thetaI=theta,
    #                 thetaC=thetaC,thetaF=thetaFinput,
    #                 addScaleParam=addScaleParam,
    #                 weightEmInf=emWeight,
    #                 weightInf=stepWeight,
    #                 verbose=verbose
    # )
    
    rownames(res$b) <- colnames(X)
    if(!missing(random)){
      rownames(res$u) <- unlist(lapply(Z, colnames))
    }
    rownames(res$monitor) <- unlist(rTermsNames)
    res$sigma <- res$monitor[,which(res$llik[1,] == max(res$llik[1,]))] # we return the ones with max llik
    res$data <- data
    res$y <- yvar
    res$partitionsX <- partitionsX
    names(res$partitions) <- unlist(rTermsNames)[1:(length(rTermsNames)-1)]
    if(!missing(random)){
      res$partitions <- lapply(res$partitions, function(x){vv= matrix(x[1,1]:x[nrow(x),ncol(x)], nrow = 1); return(vv)})
    }
    res$partitionsAll <- c(res$partitionsX,res$partitions)
    res$Dtable <- data.frame(term=names(res$partitionsAll),include=FALSE,average=FALSE,D=FALSE)
    class(res)<-c("mmec")
  }
  return(res)
}
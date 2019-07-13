mmer <- function(fixed, random, rcov, data, weights, 
                 iters=20, tolpar = 1e-03, tolparinv = 1e-06, 
                 init=NULL, constraints=NULL, method="NR", 
                 getPEV=TRUE,
                 na.method.X="exclude",
                 na.method.Y="exclude",
                 return.param=FALSE, 
                 date.warning=TRUE,
                 verbose=TRUE,reshape.output=TRUE){
  
  my.year <- 2019
  my.month <- 10 #month when the user will start to get notifications the 1st day of next month
  ### if my month = 3, user will start to get notification in april 1st (next month)
  datee <- Sys.Date()
  year.mo.day <- as.numeric(strsplit(as.character(datee),"-")[[1]])# <- as.numeric(strsplit(gsub("....-","",datee),"-")[[1]])
  your.year <- year.mo.day[1]
  your.month <- year.mo.day[2]
  ## if your month is greater than my month you are outdated
  if(date.warning){
    if(your.month > my.month & your.year >= my.year){
      # error if your month is greater and your year is smaller
      cat("Version out of date. Please update sommer to the newest version using:\ninstall.packages('sommer') in a new session\n Use the 'date.warning' argument to disable the warning message.")
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
  ## do the needed for na.method.Y and na.method.X
  dataor <- data
  provdat <- subdata(data, fixed=fixed, na.method.Y = na.method.Y,na.method.X=na.method.X)
  # print(str(provdat))
  data <- provdat$datar
  # randomization <- sample(1:nrow(data))
  # data <- data[randomization,]
  #################
  data$units <- levels(as.factor(paste("u",1:nrow(data),sep="")))
  #################
  ## get Y matrix
  response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
  responsef <- as.formula(paste(response,"~1"))
  mfna <- try(model.frame(responsef, data = data, na.action = na.pass), silent = TRUE)
  if (class(mfna) == "try-error") {
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", call. = FALSE)
  }
  mfna <- eval(mfna, data, parent.frame())
  yvar <- as.matrix(model.response(mfna))
  nt <- ncol(yvar)
  if(nt==1){colnames(yvar) <- response}
  
  #################
  ## get Zs and Ks
  if(!missing(random)){
    yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    rtermss <- apply(data.frame(yuyu),1,function(x){
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    zs <- list()
    ks <- list()
    ges <- list()
    gesI <- list()
    re_namel1 <- list()
    # bn=0.1
    bnmm <- matrix(0.1,nt,nt)+diag(.05,nt)
    # print(str(data))
    # print(str(model.matrix(~id-1,data)))
    # print(str(as(model.matrix(~id-1,data), Class = "sparseMatrix")))
    for(u in 1:length(rtermss)){
      checkvs <- grep("vs\\(",rtermss[u])
      if(length(checkvs)>0){ ## if this term is a variance structure
        # print(mm)
        ff <- eval(parse(text = rtermss[u]),data, parent.frame())# envir = data) 
        # print(nrow(ff$Z[[1]]))
        # print(length(provdat$good))
        if(nrow(ff$Z[[1]]) != length(provdat$good)){
          ## if the incidence matrix is different than the size of good very likely the user provided a matrix
          ff$Z <- lapply(ff$Z,function(xxx){as.matrix(xxx[provdat$good,])})
        }
        # print(rtermss[u])
        # print(ff$Gtc)
        re_namel1[[u]] <- ff$re_name
        zs[[u]] <- ff$Z
        ks[[u]] <- ff$K
        if(is.null(ff$Gt)){ ## initial vc if user don't provide them
          mml <- list()
          for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
            if(ff$typevc[k] == 2){div=2}else{div=1}## divisor
            # mml[[k]] <- (( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt))/div 
            mml[[k]] <- bnmm/div
          }
          ges[[u]] <- mml
        }else{ ## user provides initial values
          if(is.list(ff$Gt)){ ## user provides a list of initial values
            ges[[u]] <- ff$Gt
          }else{
            ges[[u]] <- rep(list(ff$Gt),length(ff$Z))
          }
        }
        # print(str(ff))
        # print(ff$Gtc)
        if(is.null(ff$Gtc)){ ## contraints if user don't provide them
          mml <- list()
          for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
            mm <- matrix(as.vector(ff$typevc[k]),nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
            mml[[k]] <- mm
          }
          gesI[[u]] <- mml
        }else{
          if(is.list(ff$Gtc)){# user provides a list of constraints
            gesI[[u]] <- lapply(ff$Gtc, function(x){x[lower.tri(x)] <- 0; return(x)})
          }else{
            ff$Gtc[lower.tri(ff$Gtc)] <- 0
            gesI[[u]] <- rep(list(ff$Gtc),length(ff$Z))
          }
          
        }
        # print(ff$Gtc)
        # print(gesI[[u]])
      }else{ ## if is a normal term
        zpp <- model.matrix(as.formula(paste("~",rtermss[u],"-1")), data=data)
        zs[[u]] <- list(zpp)
        nu <- ncol(zpp)
        ks[[u]] <- list(diag(1, nu,nu))
        re_namel1[[u]] <- list(rtermss[u])
        ges[[u]] <- list(bnmm) 
        mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
        gesI[[u]] <- list(mm)
        # }
      }
    }
    Z <- unlist(zs,recursive=FALSE)
    K <- unlist(ks,recursive=FALSE)
    # print(ges)
    ges <- unlist(ges,recursive=FALSE)
    gesI <- unlist(gesI,recursive=FALSE)
    re_namel1 <- unlist(re_namel1,recursive=FALSE)
  }else{re_namel1 <- character()}
  #################
  ## get Rs
  yuyur <- strsplit(as.character(rcov[2]), split = "[+]")[[1]]
  rcovtermss <- apply(data.frame(yuyur),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  rs <- list()
  gesr <- list()
  gesIr <- list()
  re_namel2 <- list()
  bnmm <- matrix(0.1,nt,nt)+diag(.05,nt)
  for(u in 1:length(rcovtermss)){
    checkvs <- grep("vs\\(",rcovtermss[u])
    if(length(checkvs)>0){ ## if this term is a variance structure
      ff <- eval(parse(text = rcovtermss[u]),data,parent.frame())# envir = data)
      rs[[u]] <- ff$Z
      re_namel2[[u]] <- ff$re_name
      if(is.null(ff$Gt)){ ## initial vc if user don't provide them
        mml <- list()
        for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
          if(ff$typevc[k] == 2){div=2}else{div=1}## divisor
          # mml[[k]] <- (( matrix(1,nt,nt) * 0 + 1) * 0.04977728 + diag(0.02488864, nt,nt))/div 
          mml[[k]] <- (bnmm/div)*5
        }
        gesr[[u]] <-  mml
      }else{gesr[[u]] <- rep(list(ff$Gt),length(ff$Z))}
      if(is.null(ff$Gtc)){ ## contraints if user don't provide them
        mml <- list()
        for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
          mm <- matrix(as.vector(ff$typevc[k]),nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
          mml[[k]] <- mm
        }
        gesIr[[u]] <- mml
      }else{
        ff$Gtc[lower.tri(ff$Gtc)] <- 0
        gesIr[[u]] <- rep(list(ff$Gtc),length(ff$Z))
      }
    }else{ ## if is a normal term
      rpp <- model.matrix(as.formula(paste("~",rcovtermss[u],"-1")), data=data)
      rs[[u]] <- list(rpp)
      re_namel2[[u]] <- rcovtermss[u]
      # gesr[[u]] <- list(( matrix(1,nt,nt) * 0 + 1) * 0.04977728 + diag(0.02488864, nt,nt))
      gesr[[u]] <- list(bnmm*5)
      mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
      gesIr[[u]] <- list(mm)
    }
  }
  R <- unlist(rs,recursive=FALSE)
  gesr <- unlist(gesr,recursive=FALSE)
  gesIr <- unlist(gesIr,recursive=FALSE)
  re_namel2 <- unlist(re_namel2,recursive=FALSE)
  
  if(!missing(random)){
    GES <- c(ges,gesr)
    GESI <- c(gesI,gesIr)
  }else{
    GES <- c(gesr)
    GESI <- c(gesIr)
  }
  
  #################
  #################
  ## get Xs
  
  yuyuf <- strsplit(as.character(fixed[3]), split = "[+-]")[[1]]
  fixedtermss <- apply(data.frame(yuyuf),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  # identify if user wants intercept or not
  test1 <- length(which(fixedtermss %in% "1"))
  test2 <- length(which(fixedtermss %in% "-1"))
  
  if(test1 == 0 & test2 == 0){ # there should be intercept
    useinter <- TRUE
  }else if(test1 > 0 & test2 == 0){# there should be intercept as well
    useinter <- TRUE
  }else{ # there's no intercept
    useinter <- FALSE
  }
  
  vsterms <- grep("vs\\(",fixedtermss)
  if(length(vsterms)>0){
    fixedvsterms <- fixedtermss[c(vsterms)]
    fixedtermss <- fixedtermss[-c(vsterms)]
    addxs <- list()
    for(o in 1:length(fixedvsterms)){
      ffx <- eval(parse(text = fixedvsterms[o]),data,parent.frame())
      addxs[[o]] <- ffx$Z[[1]]
    }
    addxs <- do.call(cbind,addxs)
  }else{fixedvsterms <- NULL}
  
  if(test1 == 0 & test2 == 0){ # there should be intercept
    fixedtermss <- c("1",fixedtermss)
  }
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  newfixed <- as.formula(paste("~",paste(fixedtermss,collapse = "+")))
  mf <- try(model.frame(newfixed, data = data, na.action = na.pass), silent = TRUE)
  mf <- eval(mf, parent.frame())
  baseX <- model.matrix(newfixed, mf)
  if(length(vsterms) > 0){baseX <- cbind(baseX,addxs)}
  qr <- qr(baseX)
  keepx <- qr$pivot[1:qr$rank]
  nameskeepx <- colnames(baseX)[keepx]
  colsdropped <- ncol(baseX) - length(nameskeepx)
  baseX <- matrix(baseX[, keepx],nrow(yvar),qr$rank)
  colnames(baseX) <- nameskeepx
  if(colsdropped > 0){
    cat(blue(paste("fixed-effect model matrix is rank deficient so dropping",colsdropped,"columns / coefficients\n")))
  }
  # print(colnames(baseX))
  xs <- list()
  gesf <- list()
  gesIf <- list()
  allfixedterms <- c(fixedtermss, fixedvsterms)
  for(u in 1:length(allfixedterms)){
    checkvs <- grep("vs\\(",allfixedterms[u])
    if(length(checkvs)>0){ ## if this term is a variance structure
      ffx <- eval(parse(text = allfixedterms[u]),data,parent.frame())# envir = data)
      findlevs <- colnames(ffx$Z[[1]])
      findlevs2 <- which(colnames(baseX) %in% findlevs)
      xpp <- as.matrix(baseX[, findlevs2])
      colnames(xpp) <- colnames(baseX)[findlevs2]
      # print(ff$Gtc)
      xs[[u]] <- list(xpp)
      if(is.null(ffx$Gt)){ ## initial vc if user don't provide them
        gesf[[u]] <- list(( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt)) 
      }else{gesf[[u]] <- rep(list(ffx$Gt),length(ffx$Z))}
      
      if(is.null(ffx$Gtc)){ ## contraints if user don't provide them
        mm <- diag(1,nt,nt); mm[lower.tri(mm)] <- 0; #mm[upper.tri(mm)] <- 2
        gesIf[[u]] <- list(mm)
      }else{gesIf[[u]] <- rep(list(ffx$Gtc),length(ffx$Z))}
    }else{ ## if is a normal term
      
      if(allfixedterms[u] == "1"){
        findlevs <- colnames(model.matrix(as.formula(paste("~",allfixedterms[u])), data=data))
      }else{
        findlevs <- colnames(model.matrix(as.formula(paste("~",allfixedterms[u],"-1")), data=data)) 
        
        check9 <- grep(":",allfixedterms[u])
        if(length(check9) > 0){
          splitby0 <- strsplit(allfixedterms[u],":")[[1]]
          splitby1 <- rep(list(splitby0), length(splitby0))
          splitby1 <- expand.grid(splitby1)
          splitby1 <- splitby1[which(apply(splitby1,1,function(x){length(unique(x))}) > 1),]
          splitby2 <- lapply(as.list(1:nrow(splitby1)), function(x){colnames(model.matrix(as.formula(paste("~", paste(as.vector(t(splitby1[x,])),collapse = ":"), "-1")), data = data))})
          matches0 <- unlist(lapply(splitby2,function(x){length(which(colnames(baseX) %in% x))}))
          best <- which(matches0 == max(matches0))[1]
          findlevs <- colnames(model.matrix(as.formula(paste("~",paste(as.vector(t(splitby1[best,])),collapse = ":"), "-1")), data = data))
        }
        check10 <- grep("[*]",allfixedterms[u])
        if(length(check10) > 0){
          splitby0 <- strsplit(allfixedterms[u],"[*]")[[1]]
          splitby1 <- rep(list(splitby0), length(splitby0))
          splitby1 <- expand.grid(splitby1)
          splitby1 <- splitby1[which(apply(splitby1,1,function(x){length(unique(x))}) > 1),]
          splitby2 <- lapply(as.list(1:nrow(splitby1)), function(x){colnames(model.matrix(as.formula(paste("~", paste(as.vector(t(splitby1[x,])),collapse = "*"), "-1")), data = data))})
          matches0 <- unlist(lapply(splitby2,function(x){length(which(colnames(baseX) %in% x))}))
          best <- which(matches0 == max(matches0))[1]
          findlevs <- colnames(model.matrix(as.formula(paste("~",paste(as.vector(t(splitby1[best,])),collapse = "*"), "-1")), data = data))
        }
        # 
      }
      # print(allfixedterms[u])
      # print(findlevs)
      findlevs2 <- which(colnames(baseX) %in% findlevs)
      xpp <- as.matrix(baseX[, findlevs2])
      colnames(xpp) <- colnames(baseX)[findlevs2]
      
      xs[[u]] <- list(xpp)
      gesf[[u]] <- list(( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt)) 
      mm <- diag(1,nt,nt); mm[lower.tri(mm)] <- 0; #mm[upper.tri(mm)] <- 2
      gesIf[[u]] <- list(mm)
    }
  }
  
  # print(str(xs))
  # print(str(gesIf))
  gesIx <- unlist(gesIf,recursive=FALSE)
  X <- unlist(xs,recursive=FALSE)
  Gx <- gesIx
  #################
  #################
  #################
  ## weights
  
  if(!missing(weights)){
    col1 <- deparse(substitute(weights))
    coco <- data[[col1]]
    ws<- coco
  }else{ws <- rep(1,nrow(yvar))}
  
  #################
  ## subset data
  if(method == "NR"){
    selected <- FALSE
  }else{selected <- TRUE}
  #################
  ## provide arguments to MNR
  if(!missing(random)){
    Z <- lapply(Z,function(x){as(x, Class = "sparseMatrix")})
  }else{Z <- list();K <- list();random=NULL}
  R <- lapply(R,function(x){as(x, Class = "sparseMatrix")})
  if(!is.null(init)){GES <- init}
  if(!is.null(constraints)){GESI <- constraints}
  re_names <- c(re_namel1,re_namel2)
  
  if(return.param){
    good <- provdat$good
    if(missing(weights)){
      args <- list(fixed=fixed, random=random, rcov=rcov, data=data, 
                        iters=iters, tolpar=tolpar, tolparinv=tolparinv, 
                        init=init, constraints=constraints, method=method, 
                        getPEV=getPEV,
                        na.method.X=na.method.X,
                        na.method.Y=na.method.Y,
                        return.param=return.param, 
                        date.warning=date.warning,
                        verbose=verbose,reshape.output=reshape.output)
    }else{
      args <- list(fixed=fixed, random=random, rcov=rcov, data=data, weights=weights, 
                        iters=iters, tolpar=tolpar, tolparinv=tolparinv, 
                        init=init, constraints=constraints, method=method, 
                        getPEV=getPEV,
                        na.method.X=na.method.X,
                        na.method.Y=na.method.Y,
                        return.param=return.param, 
                        date.warning=date.warning,
                        verbose=verbose,reshape.output=reshape.output)
    }
   
    res <- list(yvar=yvar, X=X,Gx=Gx,Z=Z,K=K,R=R,GES=GES,GESI=GESI, ws=ws,
                iters=iters, tolpar=tolpar, tolparinv=tolparinv, 
                selected=selected,getPEV=getPEV,verbose=verbose, retscaled=FALSE,
                re_names=re_names,good=good,fixedtermss=fixedtermss,args=args
    )
  }else{
    res <- .Call("_sommer_MNR",PACKAGE = "sommer",yvar, X,Gx,Z,K,R,GES,GESI, ws,
                 iters, tolpar, tolparinv,
                 selected,getPEV,verbose, FALSE)
    
    # res <- MNR(yvar, X,Gx,Z,K,R,GES,GESI, ws,
    #              iters, tolpar, tolparinv,
    #              selected,getPEV,verbose, FALSE)
    
    nslices <- dim(res$sigma)[3]
    itraits <- colnames(yvar)
    re_names <- gsub("\\(Intercept):","",re_names)
    re_names_onlyrandom <- gsub("\\(Intercept):","",re_namel1)
    dimnames(res$sigma) <- list(itraits,itraits,re_names)
    names(res$U) <- re_names_onlyrandom
    names(res$VarU) <- re_names_onlyrandom
    names(res$PevU) <- re_names_onlyrandom
    res$method <- method
    res$call <- list(fixed=fixed,random=random,rcov=rcov,
                     na.method.Y=na.method.Y,na.method.X=na.method.X)
    if(!missing(random)){
      namelist <- lapply(Z,function(x){colnames(x)})
    }else{namelist <- list()}
    # print(Gx)
    namesbeta <- lapply(as.list(1:length(X)),function(x){
      # print(Gx[[x]])
      tt <- colnames(yvar)[as.logical(apply(Gx[[x]],1,sum))]
      # print(head(X[[x]]))
      ttp <- expand.grid(tt,colnames(X[[x]]))
      # print(ttp)
      return(ttp)
    })
    # print(namesbeta)
    namesbeta <- do.call(rbind,namesbeta); 
    namelist[[length(namelist)+1]] <- namesbeta
    # res$namelist <- namelist
    if(reshape.output){
      res <- reshape_mmer(res,namelist) 
    }
    res$constraints <- GESI
    res$constraintsF <- Gx
    res$data <- data#dataor
    class(res)<-c("mmer")
  }
  
  return(res)
}

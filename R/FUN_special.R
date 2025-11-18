## small matrix constructors
unsm <- function(x, reps=NULL){
  mm <- matrix(1,x,x)
  mm[upper.tri(mm)] <- 2
  mm[lower.tri(mm)] <- 2
  if(!is.null(reps)){
    return(rep(list(mm),reps))
  }else{return(mm)}
}

unsm2 <- function (x, reps = NULL) {
  mm <- matrix(1, x, x)
  mm[upper.tri(mm)] <- 2
  mm[lower.tri(mm)] <- 0
  if (!is.null(reps)) {
    return(rep(list(mm), reps))
  }
  else {
    return(mm)
  }
}

fixm <- function(x, reps=NULL){
  mm <- matrix(3,x,x)
  if(!is.null(reps)){
    return(rep(list(mm),reps))
  }else{return(mm)}
}

covm <- function(ran1,ran2, thetaC=NULL, theta=NULL){
  if( ncol(ran1$Z[[1]]) != ncol(ran2$Z[[1]]) ){stop("Matrices of the two random effects should have the same dimensions",call. = FALSE)}
  ran1$Z[[2]] <- ran2$Z[[1]]
  if(is.null(thetaC)){
    ran1$thetaC <- unsm(2); 
  }else{ran1$thetaC <- thetaC}
  ran1$thetaC[lower.tri(ran1$thetaC)] = 0 # lower.tri must be 0
  colnames(ran1$thetaC) <- rownames(ran1$thetaC) <- c("ran1","ran2")
  if(is.null(theta)){
    ran1$theta <- diag(2) * 0.15 + matrix(0.015, 2, 2)
  }else{ran1$theta <- theta}
  nvc <- length(which(thetaC > 0))
  ran1$thetaF <- diag(nvc) # n x n, where n is number of vc to estimate
  ran1$sp <- rep(0,nvc) # rep 0 n times, where n is number of vc to estimate
  return(ran1)
}

replace.values <- function(Values,Search,Replace){
  dd0 <- data.frame(Values)
  vv <- which(Values%in%Search)
  dd <- data.frame(Search,Replace)
  rownames(dd) <- Search
  dd0[vv,"Values"] <- as.character(dd[Values[vv],"Replace"])
  return(dd0[,1])
}

myformula <- function(x){
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  yuyuf <- strsplit(as.character(x[3]), split = "[+]")[[1]]
  termss <- apply(data.frame(yuyuf),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  newtermss <- apply(data.frame(yuyuf),1,function(y){
    newy <- expi(y)
    if(length(newy) > 0){
      newy <- gsub(",.*","",newy)
    }else{newy <- y}
    return(newy)
  })
  resp <- strsplit(as.character(x[2]), split = "[+]")[[1]]
  newx <- paste(resp, "~",paste(newtermss,collapse = "+"))
  return(newx)
}

##############
## na.methods

subdata <- function(data,fixed,na.method.Y=NULL,na.method.X=NULL){
  
  # silently change all columns that are defined as character into factors
  # columnTypes <- unlist(lapply(data, class))
  # columnTypesC <- which(columnTypes == "character")
  # if(length(columnTypesC) > 0){ # if there's character types change them to factor
  #   for(cti in 1:length(columnTypesC)){
  #     data[,cti] <- as.factor(data[,cti])
  #   }
  # }
  ####
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
  responsef <- as.formula(paste(response,"~1"))
  mfna <- try(model.frame(responsef, data = data, na.action = na.pass), silent = TRUE)
  if (is(mfna, "try-error") ) { # class(mfna) == "try-error"
    stop("Please provide the 'data' argument for your specified variables.\nYou may be specifying some variables in your model not present in your dataset.", call. = FALSE)
  }
  mfna <- eval(mfna, parent.frame())
  yvar <- as.matrix(model.response(mfna))
  nt <- ncol(yvar)
  good <- 1:nrow(data)
  if(nt==1){colnames(yvar) <- response}
  if(na.method.Y=="include"){
    touse <- colnames(yvar)
    for(i in 1:length(touse)){
      use <- touse[i]
      # print(iname)
      data[,use] <- imputev(data[,use])
    }
  }else if(na.method.Y=="include2"){
    tlist <- list()
    touse <- colnames(yvar)
    for(i in 1:length(touse)){
      # print(touse[i])
      use <- touse[i]
      vivi <- as.vector(data[,use])
      tlist[[i]] <- which(!is.na(vivi))
    }
    # print(tlist)
    good <- sort(unique(unlist(tlist)))
    data <- data[good,]
    for(i in 1:length(touse)){
      use <- touse[i]
      # print(iname)
      data[,use] <- imputev(data[,use])
    }
  }else if(na.method.Y=="exclude"){
    tlist <- list()
    touse <- colnames(yvar)
    for(i in 1:length(touse)){
      # print(touse[i])
      use <- touse[i]
      vivi <- as.vector(data[,use])
      tlist[[i]] <- which(!is.na(vivi))
    }
    # print(tlist)
    if(length(tlist)==1){ #only one trait
      good <- tlist[[1]]
    }else{#more than one trait
      good <- Reduce(intersect,tlist)
    }
    data <- data[good,]
  }else{stop("na.method.Y not recognized")}
  data <- data.frame(data)
  
  ##########
  ## na.method x
  yuyu <- strsplit(as.character(fixed[3]), split = "[+]")[[1]]
  xtermss <- apply(data.frame(yuyu),1,function(x){
    strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
  })
  xtermss2 <- apply(data.frame(xtermss),1,function(x){gsub(",.*","",expi2(x))})
  xtermss2[which(xtermss2 == "")] <- xtermss[which(xtermss2 == "")]
  
  xtermss2 <- intersect(colnames(data),xtermss2) # only focus on the terms that are in teh dataset so we can skip overlay and weird vs structures
  
  # print(xtermss2)
  if(length(xtermss2) > 0){
    mycl <- as.vector(unlist(lapply(data.frame(data[,xtermss2]),class)))
    
    if(na.method.X=="include"){
      touse <- xtermss2
      for(i in 1:length(touse)){
        use <- touse[i]
        usecl <- mycl[i]
        if(usecl == "factor"){data[,use] <- as.factor(imputev(data[,use]))}else{data[,use] <- imputev(data[,use])}
      }
    }else if(na.method.X=="exclude"){
      tlist <- list()
      touse <- xtermss2
      for(i in 1:length(touse)){
        # print(touse[i])
        use <- touse[i]
        vivi <- as.vector(data[,use])
        tlist[[i]] <- which(!is.na(vivi))
      }
      # print(tlist)
      if(length(tlist)==1){ #only one trait
        good <- tlist[[1]]
      }else{#more than one trait
        good <- Reduce(intersect,tlist)
      }
      data <- data[good,]
    }else{stop("na.method.Y not recognized")}
    data <- data.frame(data)
  }
  
  return(list(datar=data,good=good))
  
}

###############
## VS structures for mmec

H <- function(timevar=NULL, idvar=NULL, response=NULL, Gu=NULL){
  if(is.null(timevar) ){stop("Please provide the timevar argument.", call. = FALSE)}
  if(is.null(idvar) ){stop("Please provide the idvar argument.", call. = FALSE)}
  if(is.null(response) ){stop("Please provide the response argument.", call. = FALSE)}
  
  dtx <- data.frame(timevar=timevar, idvar=idvar, v.names=response)
  dtx2 <- aggregate(v.names~timevar+idvar, data=dtx, FUN=mean, na.rm=TRUE)
  wide <- reshape(dtx2, direction = "wide", idvar = "idvar",
                  timevar = "timevar", v.names = "v.names", sep= "_")
  rowNamesWide <-  wide[,1]
  rownames(wide) <- rowNamesWide
  wide <- wide[,-1]
  # if user doesn't provide the a Gu we impute simply and use the correlation matrix as a Gu
  if(is.null(Gu)){ 
    X <- apply(wide, 2, imputev)
    Gu <- cor(t(X))
  }else{
    Gu = cov2cor(Gu)
  } 
  # impute missing data using a relationship matrix 
  if(is.null(rownames(Gu))){stop("Gu needs to have row names.", call. = FALSE)}
  if(is.null(colnames(Gu))){stop("Gu needs to have column names.", call. = FALSE)}
  for(iEnv in 1:ncol(wide)){ # iEnv=1
    withData <- which(!is.na(wide[,iEnv]))
    withoutData <- which(is.na(wide[,iEnv]))
    imputationVector <- as.numeric(Gu[as.character(rowNamesWide),as.character(rowNamesWide[withData])] %*% as.matrix(wide[withData,iEnv]))
    wide[,iEnv] <- imputationVector  # wide[withoutData,iEnv] <- imputationVector[withoutData]
    # scaleFactor=imputationVector[withData[1]] / wide[withData[1],iEnv]
  }
  colnames(wide) <- gsub("v.names_","", colnames(wide))
  return(wide)
}

atm <- function(x, levs, thetaC=NULL, theta=NULL){
  if(is.matrix(x)){
    dummy <- x
    dummy <- dummy[,levs]
    m0 <- rep(0,ncol(dummy))
    names(m0) <- levels(as.factor(colnames(dummy)))#as.character(unique(x))
    if(missing(levs)){levs <- names(m0)}
    m0[levs] <- 1
    mm <- Diagonal(m0)
    colnames(mm) <- rownames(mm) <- colnames(dummy)
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- Matrix(x,ncol=1); colnames(dummy) <- namess
      dummy <- dummy[,levs]
      m0 <- rep(0,ncol(dummy))
      names(m0) <- levels(as.factor(colnames(dummy)))#as.character(unique(x))
      if(missing(levs)){levs <- names(m0)}
      m0[levs] <- 1
      mm <- Diagonal(m0)
      colnames(mm) <- rownames(mm) <- colnames(dummy)
    }else{
      dummy <- x
      dummy <- Matrix::sparse.model.matrix(~dummy-1,na.action = na.pass)
      colnames(dummy) <- gsub("dummy","",colnames(dummy))
      dummy <- dummy[,levs, drop=FALSE]
      m0 <- rep(0,ncol(dummy))
      names(m0) <- levels(as.factor(colnames(dummy)))#as.character(unique(x))
      if(missing(levs)){levs <- names(m0)}
      m0[levs] <- 1
      mm <- diag(m0)
      colnames(mm) <- rownames(mm) <- colnames(dummy)
    }
  }
  bnmm <- mm*0.15
  if(nrow(bnmm) > 1){
    bnmm[upper.tri(bnmm)]=bnmm[upper.tri(bnmm)]/2
  }
  if(!is.null(thetaC)){
    mm <- thetaC
    colnames(mm) <- rownames(mm) <- colnames(dummy)
  }
  if(!is.null(theta)){
    bnmm <- theta
    colnames(bnmm) <- rownames(bnmm) <- colnames(dummy)
  }
  mm[lower.tri(mm)]=0
  return(list(Z=dummy,thetaC=mm, theta=bnmm))
}
csm <- function(x,mm, thetaC=NULL, theta=NULL){
  if(is.matrix(x)){
    mm <- mm
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- Matrix(x,ncol=1); colnames(dummy) <- namess
    }else{
      dummy <- x
      levs <- na.omit(unique(dummy))
      if(length(levs) > 1){
        dummy  <- Matrix::sparse.model.matrix(~dummy-1,na.action = na.pass)
        colnames(dummy) <- gsub("dummy","",colnames(dummy))
      }else{
        vv <- which(!is.na(dummy));
        dummy <- matrix(0,nrow=length(dummy))
        dummy[vv,] <- 1; colnames(dummy) <- levs
      }
    }
    mm <- mm
  }
  # mm[lower.tri(mm)] <- 0
  colnames(mm) <- rownames(mm) <- colnames(dummy)
  bnmm <- mm*0.15
  if(nrow(bnmm) > 1){
    bnmm[upper.tri(bnmm)]=bnmm[upper.tri(bnmm)]/2
  }
  if(!is.null(thetaC)){
    mm <- thetaC
    colnames(mm) <- rownames(mm) <- colnames(dummy)
  }
  if(!is.null(theta)){
    bnmm <- theta
    colnames(bnmm) <- rownames(bnmm) <- colnames(dummy)
  }
  mm[lower.tri(mm)]=0
  return(list(Z=dummy,thetaC=mm,theta=bnmm))
}
dsm <- function(x, thetaC=NULL, theta=NULL){
  if(is.matrix(x)){
    dummy <- x
    mm <- diag(1,ncol(x))
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <-  as(as(as( Matrix(x,ncol=1) ,  "dMatrix"), "generalMatrix"), "CsparseMatrix") # as(Matrix(x,ncol=1), Class = "dgCMatrix"); 
      colnames(dummy) <- namess
      mm <- diag(ncol(dummy));
    }else{
      dummy <- x
      levs <- na.omit(unique(dummy))
      if(length(levs) > 1){
        dummy  <- Matrix::sparse.model.matrix(~dummy-1,na.action = na.pass)
        colnames(dummy) <- gsub("dummy","",colnames(dummy))
      }else{
        vv <- which(!is.na(dummy));
        dummy <- matrix(0,nrow=length(dummy))
        dummy[vv,] <- 1; colnames(dummy) <- levs
      }
      mm <- diag(1,ncol(dummy))
    }
  }
  colnames(mm) <- rownames(mm) <- colnames(dummy)
  bnmm <- mm*0.15
  if(nrow(bnmm) > 1){
    bnmm[upper.tri(bnmm)]=bnmm[upper.tri(bnmm)]/2
  }
  if(!is.null(thetaC)){
    mm <- thetaC
  }
  if(!is.null(theta)){
    bnmm <- theta
  }
  colnames(mm) <- rownames(mm) <- colnames(dummy)
  colnames(bnmm) <- rownames(bnmm) <- colnames(dummy)
  mm[lower.tri(mm)]=0
  return(list(Z=dummy,thetaC=mm, theta=bnmm))
}
usm <- function(x, thetaC=NULL, theta=NULL){
  # namx <- as.character(substitute(list(x)))[-1L]
  if(is.matrix(x)){
    dummy <- x
    mm <- unsm(ncol(dummy))
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- Matrix(x,ncol=1); colnames(dummy) <- namess
    }else{
      dummy <- x
      levs <- na.omit(unique(dummy))
      if(length(levs) > 1){
        dummy  <- Matrix::sparse.model.matrix(~dummy-1,na.action = na.pass)
        colnames(dummy) <- gsub("dummy","",colnames(dummy))
      }else{
        vv <- which(!is.na(dummy));
        dummy <- Matrix(0,nrow=length(dummy))
        dummy[vv,] <- 1; colnames(dummy) <- levs
      }
    }
    mm <- unsm(ncol(dummy))
  }
  
  colnames(mm) <- rownames(mm) <- colnames(dummy)
  bnmm <- diag(ncol(mm))*.05 + matrix(.1,ncol(mm),ncol(mm))
  if(!is.null(thetaC)){
    mm <- thetaC
  }
  if(!is.null(theta)){
    bnmm <- theta
  }
  colnames(mm) <- rownames(mm) <- colnames(dummy)
  colnames(bnmm) <- rownames(bnmm) <- colnames(dummy)
  mm[lower.tri(mm)]=0
  return(list(Z=dummy,thetaC=mm,theta=bnmm))
}
ism <- function(x, thetaC=NULL, theta=NULL){
  if(class(x)[1] %in% c("dgCMatrix","matrix") ){
    dummy <-  as(as(as( x ,  "dMatrix"), "generalMatrix"), "CsparseMatrix") # as(x, Class="dgCMatrix")
    mm <- diag(1)#,ncol(x))
  }else{ # if user provides a vector
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- Matrix(x,ncol=1); colnames(dummy) <- namess
      mm <- diag(1);
    }else{
      dummy <- x
      levs <- na.omit(unique(dummy))
      if(length(levs) > 1){
        dummy  <- Matrix::sparse.model.matrix(~dummy-1,na.action = na.pass)
        if(!is(class(dummy), "dgCMatrix")){
          dummy <-  as(as(as( dummy ,  "dMatrix"), "generalMatrix"), "CsparseMatrix") # as(dummy, Class="dgCMatrix")
        }
        colnames(dummy) <- gsub("dummy","",colnames(dummy))
      }else{
        vv <- which(!is.na(dummy));
        dummy <- Matrix(0,nrow=length(dummy))
        dummy[vv,] <- 1; colnames(dummy) <- levs
      }
      mm <- diag(1)
    }
  }
  colnames(mm) <- rownames(mm) <- "mu"# colnames(dummy)
  bnmm <- mm*.15
  if(nrow(bnmm) > 1){
    bnmm[upper.tri(bnmm)]=bnmm[upper.tri(bnmm)]/2
  }
  if(!is.null(thetaC)){
    mm <- thetaC
  }
  if(!is.null(theta)){
    bnmm <- theta
  }
  mm[lower.tri(mm)]=0
  return(list(Z=dummy,thetaC=mm, theta=bnmm))
}




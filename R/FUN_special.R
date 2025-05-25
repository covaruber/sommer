## small matrix constructors
unsm <- function(x, reps=NULL){
  mm <- matrix(1,x,x)
  mm[upper.tri(mm)] <- 2
  mm[lower.tri(mm)] <- 2
  if(!is.null(reps)){
    return(rep(list(mm),reps))
  }else{return(mm)}
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

stackTrait <- function (data, traits) 
{
  dataScaled <- data
  traits <- intersect(traits, colnames(data))
  idvars <- setdiff(colnames(data), traits)
  for (iTrait in traits) {
    dataScaled[, iTrait] <- scale(dataScaled[, iTrait])
  }
  columnTypes <- unlist(lapply(data[idvars], class))
  columnTypes <- columnTypes[which(columnTypes %in% c("factor", 
                                                      "character"))]
  idvars <- intersect(idvars, names(columnTypes))
  data2 <- reshape(data[, c(idvars, traits)], idvar = idvars, 
                   varying = traits, timevar = "trait", times = traits, 
                   v.names = "value", direction = "long")
  data2Scaled <- reshape(dataScaled[, c(idvars, traits)], idvar = idvars, 
                         varying = traits, timevar = "trait", times = traits, 
                         v.names = "value", direction = "long")
  data2 <- as.data.frame(data2)
  data2$valueS <- as.vector(unlist(data2Scaled$value))
  rownames(data2) <- NULL
  varG <- cov(data[, traits], use = "pairwise.complete.obs")
  mu <- apply(data[, traits], 2, mean, na.rm = TRUE)
  # x
  return(list(long = data2, varG = varG, mu = mu))
}

add.diallel.vars <- function(df, par1="Par1", par2="Par2",sep.cross="-"){
  # Dummy variables for selfs, crosses, combinations
  df[,"is.cross"] <- ifelse(df[,par1] == df[,par2], 0, 1)
  df[,"is.self"] <- ifelse(df[,par1] == df[,par2], 1, 0)
  df[,"cross.type"] <- ifelse(as.character(df[,par1]) < as.character(df[,par2]), -1,
                              ifelse(as.character(df[,par1]) == as.character(df[,par2]), 0, 1))
  # Dummy variable for the combinations, ignoring the reciprocals
  df[,"cross.id"]<-factor(ifelse(as.character(df[,par1]) <= as.character(df[,par2]),
                                 paste(df[,par1], df[,par2], sep =sep.cross),
                                 paste(df[,par2], df[,par1], sep =sep.cross)) )
  return(df)
}

overlay<- function (..., rlist = NULL, prefix = NULL, sparse=FALSE){
  init <- list(...) # init <- list(DT$femalef,DT$malef)
  ## keep track of factor variables
  myTypes <- unlist(lapply(init,class))
  init0 <- init
  ##
  init <- lapply(init, as.character)
  namesInit <- as.character(substitute(list(...)))[-1L] # names <- c("femalef","malef")
  dat <- as.data.frame(do.call(cbind, init))
  dat <- as.data.frame(dat)
  ## bring back the levels
  for(j in 1:length(myTypes)){
    if(myTypes[j]=="factor"){
      levels(dat[,j]) <- c(levels(dat[,j]),setdiff(levels(init0[[j]]),levels(dat[,j]) ))
    }
  }
  ##
  if (is.null(dim(dat))) {
    stop("Please provide a data frame to the overlay function, not a vector.\\n",
         call. = FALSE)
  }
  if (is.null(rlist)) {
    rlist <- as.list(rep(1, dim(dat)[2]))
  }
  ss1 <- colnames(dat)
  dat2 <- as.data.frame(dat[, ss1])
  head(dat2)
  colnames(dat2) <- ss1
  femlist <- list()
  S1list <- list()
  for (i in 1:length(ss1)) {
    femlist[[i]] <- ss1[i]
    dat2[, femlist[[i]]] <- as.factor(dat2[, femlist[[i]]])
    if(sparse){
      S1 <- Matrix::sparse.model.matrix(as.formula(paste("~", femlist[[i]],
                                                         "-1")), dat2)
    }else{
      S1 <- model.matrix(as.formula(paste("~", femlist[[i]],
                                          "-1")), dat2)
    }
    colnames(S1) <- gsub(femlist[[i]], "", colnames(S1))
    S1list[[i]] <- S1
  }
  levo <- sort(unique(unlist(lapply(S1list, function(x) {
    colnames(x)
  }))))
  if(sparse){
    S3 <- Matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
  }else{
    S3 <- matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
  }
  
  rownames(S3) <- rownames(dat2)
  colnames(S3) <- levo
  for (i in 1:length(S1list)) {
    if (i == 1) {
      S3[rownames(S1list[[i]]), colnames(S1list[[i]])] <- S1list[[i]] *
        rlist[[i]]
    }
    else {
      S3[rownames(S1list[[i]]), colnames(S1list[[i]])] <- S3[rownames(S1list[[i]]),
                                                             colnames(S1list[[i]])] + (S1list[[i]][rownames(S1list[[i]]),
                                                                                                   colnames(S1list[[i]])] * rlist[[i]])
    }
  }
  if (!is.null(prefix)) {
    colnames(S3) <- paste(prefix, colnames(S3), sep = "")
  }
  attr(S3,"variables") <- namesInit
  return(S3)
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
redmm <- function (x, M = NULL, Lam=NULL, nPC=50, cholD=FALSE, returnLam=FALSE) {
  
  if(system.file(package = "RSpectra") == ""){
    stop("Please install the RSpectra package to use the redmm() function.",call. = FALSE)
  }else{
    requireNamespace("RSpectra",quietly=TRUE)
  }
  
  if(is.null(M)){
    # stop("M cannot be NULL. We need a matrix of features that defines the levels of x")
    smd <- RSpectra::svds(x, k=nPC, which = "LM")
    if(is.null(Lam)){
      Lam0 <- smd$u
      Lam = Lam0[,1:min(c(nPC,ncol(x))), drop=FALSE]
      rownames(Lam) <- rownames(x)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }else{
      Lam0=Lam
      Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
      rownames(Lam) <- rownames(M)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }
    Zstar <- Lam
  }else{
    
    if (inherits(x, "dgCMatrix") | inherits(x, "matrix")) {
      notPresentInM <- setdiff(colnames(Z),rownames(M))
      notPresentInZ <- setdiff(rownames(M),colnames(x))
    }else{
      notPresentInM <- setdiff(unique(x),rownames(M))
      notPresentInZ <- setdiff(rownames(M),unique(x))
    }
    if(is.null(Lam)){ # user didn't provide a Lambda matrix
      if(nPC == 0){ # user wants to use the full marker matrix
        Lam <- Lam0 <- M
      }else{ # user wants to use the PCA method
        nPC <- min(c(nPC, ncol(M)))
        if(cholD){
          smd <- try(chol(M) , silent = TRUE)
          if(inherits(smd, "try-error")){smd <- try(chol((M+diag(1e-5,nrow(M),nrow(M))) ) , silent = TRUE)}
          Lam0 = t(smd)
        }else{
          smd <- RSpectra::svds(M, k=nPC, which = "LM")
          Lam0 <- smd$u
        }
        Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
        rownames(Lam) <- rownames(M)
        colnames(Lam) <- paste0("nPC",1:nPC)
      }
    }else{ # user provided it's own Lambda matrix
      Lam0=Lam
      Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
      rownames(Lam) <- rownames(M)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }
  }
  if (inherits(x, "dgCMatrix") | inherits(x, "matrix")) {
    Z <- x
  }else{
    if (!is.character(x) & !is.factor(x)) {
      namess <- as.character(substitute(list(x)))[-1L]
      Z <- Matrix(x, ncol = 1)
      colnames(Z) <- namess
    }else {
      dummy <- x
      levs <- na.omit(unique(dummy))
      if (length(levs) > 1) {
        Z <- Matrix::sparse.model.matrix(~dummy - 1, na.action = na.pass)
        colnames(Z) <- gsub("dummy", "", colnames(Z))
      } else {
        vv <- which(!is.na(dummy))
        Z <- Matrix(0, nrow = length(dummy))
        Z[vv, ] <- 1
        colnames(Z) <- levs
      }
    }
  }
  
  if(is.null(M)){
    Zstar <- Lam
  }else{
    Zstar <- as.matrix(Z %*% Lam[colnames(Z),])
  }
  
  if(returnLam){
    return(list(Z = Zstar, Lam=Lam, Lam0=Lam0)) 
  }else{return(Zstar)}
  
}

rrm <- function (x = NULL, H = NULL, nPC = 2, returnGamma = FALSE, cholD = TRUE) 
{
  if (is.null(x)) {
    stop("Please provide the x argument.", call. = FALSE)
  }
  if (is.null(H)) {
    stop("Please provide the x argument.", call. = FALSE)
  }
  Y <- apply(H, 2, imputev)
  Sigma <- cov(scale(Y, scale = TRUE, center = TRUE))
  Sigma <- as.matrix(nearPD(Sigma)$mat)
  if (cholD) {
    Gamma <- t(chol(Sigma))
  }
  else {
    U <- svd(Sigma)$u
    D <- diag(svd(Sigma)$d)
    Gamma <- U %*% sqrt(D)
    rownames(Gamma) <- colnames(Gamma) <- rownames(Sigma)
  }
  colnamesGamma <- colnames(Gamma)
  rownamesGamma <- rownames(Gamma)
  Gamma <- Gamma[, 1:nPC, drop = FALSE]
  colnames(Gamma) <- colnamesGamma[1:nPC]
  rownames(Gamma) <- rownamesGamma
  rownames(Gamma) <- gsub("v.names_", "", rownames(Gamma))
  colnames(Gamma) <- paste("PC", 1:ncol(Gamma), sep = "")
  dtx <- data.frame(timevar = x)
  dtx$index <- 1:nrow(dtx)
  Z <- Matrix::sparse.model.matrix(~timevar - 1, na.action = na.pass, 
                                   data = dtx)
  colnames(Z) <- gsub("timevar", "", colnames(Z))
  Zstar <- Z %*% Gamma[colnames(Z), ]
  Zstar <- as.matrix(Zstar)
  rownames(Z) <- NULL
  if (returnGamma) {
    return(list(Gamma = Gamma, H = H, Sigma = Sigma, Zstar = Zstar))
  }
  else {
    return(Zstar)
  }
}

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
    X <- apply(wide, 2, sommer::imputev)
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




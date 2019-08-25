list2usmat <- function(sigmaL){
  
  f <- function(n, x){
    res <- ((n*(n-1))/2 + n) - x
    if(res < 0){res <- 100}
    return(res)
  }
  if(is.list(sigmaL)){
    ss <- unlist(sigmaL)
  }else{ss <- sigmaL}
  
  x <- length(ss)
  n <- round(optimize(f, c(1, 50), tol = 0.0001, x=x)$minimum)
  mss <- matrix(NA,n,n)
  mss[upper.tri(mss,diag = TRUE)] <- ss
  mss[lower.tri(mss)] <- t(mss[upper.tri(mss)])
  return(mss)
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

reshape_mmer <- function(object, namelist){
  
  nt <- nrow(as.matrix(object$sigma[,,1]))
  ntn <- attr(object$sigma, "dimnames")[[2]]
  nre <- length(object$U)
  nren <- attr(object$sigma, "dimnames")[[3]]
  U2 <- list()#vector("list",nre)
  VarU2 <- list()
  PevU2 <- list()
  
  if(nre > 0){
    for(ire in 1:nre){
      
      N <- nrow(object$U[[ire]])
      utlist <- list()# vector("list",nt)
      varutlist <- list()# vector("list",nt)
      pevutlist <- list()# vector("list",nt)
      
      # pick <- numeric()
      # for(it in 1:nt){
      #   pick <- c(pick,seq(it,N,nt))
      # }
      provus <- matrix(object$U[[ire]],ncol=nt,byrow = T)
      
      for(it in 1:nt){
        pick <- seq(it,N,nt)
        pit <- ntn[it]
        utlist[[pit]] <- provus[,it] # object$U[[ire]][pick,]
        names(utlist[[ pit ]]) <- namelist[[ire]]
        # print(length(object$VarU[[ire]]))
        if(length(object$VarU[[ire]]) > 0){
          varutlist[[ pit ]] <- as.matrix(object$VarU[[ire]][pick,pick])
          rownames(varutlist[[ pit ]]) <- colnames(varutlist[[ pit ]]) <- namelist[[ire]]
        }else{varutlist[[ pit ]] <- varutlist[[ pit ]] }
        if(length(object$PevU[[ire]]) > 0){
          pevutlist[[ pit ]] <-  as.matrix(object$PevU[[ire]][pick,pick])
          rownames(pevutlist[[ pit ]]) <- colnames(pevutlist[[ pit ]]) <- namelist[[ire]]
        }else{pevutlist[[ pit ]] <-pevutlist[[ pit ]]}
      }
      
      U2[[ nren[ire] ]] <- utlist
      VarU2[[ nren[ire] ]] <- varutlist
      PevU2[[ nren[ire] ]] <- pevutlist
    }
    
    object$U <- U2; U2<-NULL
    object$VarU <- VarU2; VarU2<-NULL
    object$PevU <- PevU2; PevU2<-NULL
  }
  
  
  
  N <- nrow(object$Vi)
  pick <- numeric()
  for(it in 1:nt){
    pick <- c(pick,seq(it,N,nt))
  }
  
  object$Vi <- object$Vi[pick,pick]
  
  # object$constraintsF
  # object$Beta <- matrix(object$Beta,ncol=nt,byrow = F)
  # print(namelist[[length(namelist)]])
  # print(object$Beta)
  # rownames(object$Beta) <- namelist[[length(namelist)]]
  object$Beta <- cbind(namelist[[length(namelist)]],object$Beta)
  colnames(object$Beta) <- c("Trait","Effect","Estimate")
  
  object$fitted <- matrix(object$fitted,ncol=nt, byrow=TRUE)
  object$residuals <- matrix(object$residuals,ncol=nt, byrow = TRUE)
  
  MyArray <- object$sigma
  object$sigma <- lapply(seq(dim(MyArray)[3]), function(x) MyArray[ , , x])
  object$sigma <- lapply(object$sigma,function(x){mm <- as.matrix(x);colnames(mm)<-rownames(mm) <- ntn;return(mm)})
  names(object$sigma) <- nren
  
  MyArray <- object$sigma_scaled
  object$sigma_scaled <- lapply(seq(dim(MyArray)[3]), function(x) MyArray[ , , x])
  object$sigma_scaled <- lapply(object$sigma_scaled,function(x){mm <- as.matrix(x);colnames(mm)<-rownames(mm) <- ntn;return(mm)})
  names(object$sigma_scaled) <- nren
  # colnames(object$Beta) <- ntn
  return(object)
}

overlay<- function (..., rlist = NULL, prefix = NULL){
  init <- list(...)
  init <- lapply(init, as.character)
  names <- as.character(substitute(list(...)))[-1L]
  dat <- as.data.frame(do.call(cbind, init))
  dat <- as.data.frame(dat)
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
    S1 <- model.matrix(as.formula(paste("~", femlist[[i]], 
                                        "-1")), dat2)
    colnames(S1) <- gsub(femlist[[i]], "", colnames(S1))
    S1list[[i]] <- S1
  }
  levo <- sort(unique(unlist(lapply(S1list, function(x) {
    colnames(x)
  }))))
  S3 <- matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
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
  return(S3)
}

##############
## na.methods

subdata <- function(data,fixed,na.method.Y=NULL,na.method.X=NULL){
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  response <- strsplit(as.character(fixed[2]), split = "[+]")[[1]]
  responsef <- as.formula(paste(response,"~1"))
  mfna <- try(model.frame(responsef, data = data, na.action = na.pass), silent = TRUE)
  if (class(mfna) == "try-error") {
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


############## 
## VS structure
at <- function(x, levs){
  if(is.matrix(x)){
    dummy <- x
    m0 <- rep(0,ncol(dummy))
    names(m0) <- levels(as.factor(colnames(dummy)))#as.character(unique(x))
    if(missing(levs)){levs <- names(m0)}
    m0[levs] <- 1
    mm <- diag(m0)
    colnames(mm) <- rownames(mm) <- colnames(dummy)
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- matrix(x,ncol=1); colnames(dummy) <- namess
      m0 <- rep(0,ncol(dummy))
      names(m0) <- levels(as.factor(colnames(dummy)))#as.character(unique(x))
      if(missing(levs)){levs <- names(m0)}
      m0[levs] <- 1
      mm <- diag(m0)
      colnames(mm) <- rownames(mm) <- colnames(dummy)
    }else{
      dummy <- x
      dummy <- model.matrix(~dummy-1,na.action = na.pass)
      colnames(dummy) <- gsub("dummy","",colnames(dummy))
      m0 <- rep(0,ncol(dummy))
      names(m0) <- levels(as.factor(colnames(dummy)))#as.character(unique(x))
      if(missing(levs)){levs <- names(m0)}
      m0[levs] <- 1
      mm <- diag(m0)
      colnames(mm) <- rownames(mm) <- colnames(dummy)
    }
  }
  return(list(dummy,mm))
}
cs <- function(x,mm){
  if(is.matrix(x)){
    mm <- mm
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- matrix(x,ncol=1); colnames(dummy) <- namess
    }else{
      dummy <- x
      levs <- na.omit(unique(dummy))
      if(length(levs) > 1){
        dummy  <- model.matrix(~dummy-1,na.action = na.pass)
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
  return(list(dummy,mm))
}
ds <- function(x){
  if(is.matrix(x)){
    dummy <- x
    mm <- diag(1,ncol(x))
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- matrix(x,ncol=1); colnames(dummy) <- namess
      mm <- diag(ncol(dummy)); 
    }else{
      dummy <- x
      levs <- na.omit(unique(dummy))
      if(length(levs) > 1){
        dummy  <- model.matrix(~dummy-1,na.action = na.pass)
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
  return(list(dummy,mm))
}
us <- function(x){
  # namx <- as.character(substitute(list(x)))[-1L]
  if(is.matrix(x)){
    dummy <- x
    mm <- unsm(ncol(dummy))
  }else{
    if(!is.character(x) & !is.factor(x)){
      namess <- as.character(substitute(list(x)))[-1L]
      dummy <- matrix(x,ncol=1); colnames(dummy) <- namess
    }else{
      dummy <- x
      levs <- na.omit(unique(dummy))
      if(length(levs) > 1){
        dummy  <- model.matrix(~dummy-1,na.action = na.pass)
        colnames(dummy) <- gsub("dummy","",colnames(dummy))
      }else{
        vv <- which(!is.na(dummy)); 
        dummy <- matrix(0,nrow=length(dummy))
        dummy[vv,] <- 1; colnames(dummy) <- levs
      }
    }
    mm <- unsm(ncol(dummy))
  }
  colnames(mm) <- rownames(mm) <- colnames(dummy) 
  return(list(dummy,mm))
}
unsm <- function(x, reps=NULL){
  mm <- matrix(1,x,x)
  mm[upper.tri(mm)] <- 2
  mm[lower.tri(mm)] <- 2
  if(!is.null(reps)){
    return(rep(list(mm),reps))
  }else{return(mm)}
}
uncm <- function(x, reps=NULL){
  mm <- matrix(2,x,x)
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
fcm <- function(x, reps=NULL){
  mm <- diag(x)
  mm <- mm[,which(apply(mm,2,sum) > 0)]
  mm <- as.matrix(mm)
  if(!is.null(reps)){
    return(rep(list(mm),reps))
  }else{return(mm)}
}

vs <- function(..., Gu=NULL, Gt=NULL, Gtc=NULL){
  
  ## ... list of structures to define the random effect
  ## Gu the known covariance matrix of the vs
  ## Gt the multitrait structure and constraints for it
  ## Gtc the initial values for the var-cov components
  init <- list(...)
  namess <- as.character(substitute(list(...)))[-1L]
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  
  namess2 <- apply(data.frame(namess),1,function(x){
    newx <- expi(x); if(length(newx)==0){newx<-""}
    newx <- gsub(",.*","",newx)
    return(newx)
  })
  namess2[which(namess2 == "")] <- namess[which(namess2 == "")]
  ref_name <- namess2[length(namess2)]
  if("units" %in% namess2){
    is.residual =TRUE
  }else{is.residual=FALSE}
  ### get the data
  init2 <- list()
  for(i in 1:length(init)){
    if(is.list(init[[i]])){ ## if it comes from a ds, us, cs function
      init2[[i]] <- init[[i]]
    }else{ # is a single vector with numbers or characters, ...
      if(is.matrix(init[[i]])){
        mm=diag(ncol(init[[i]])); rownames(mm) <- colnames(mm) <- colnames(init[[i]])
        init2[[i]] <- list(x=init[[i]],mm)
      }else{
        dummy <- init[[i]]
        if(!is.character(dummy) & !is.factor(dummy)){
          dummy <- matrix(dummy,ncol=1)
          colnames(dummy) <- namess2[i]
          mm=diag(1); rownames(mm) <- colnames(mm) <- namess2[i]
        }else{
          levs <- na.omit(unique(dummy))
          if(length(levs) > 1){
            dummy  <- model.matrix(~dummy-1,na.action = na.pass)
          }else{
            vv <- which(!is.na(dummy)); 
            dummy <- matrix(0,nrow=length(dummy))
            dummy[vv,] <- 1; colnames(dummy) <- levs
          }
          colnames(dummy) <- gsub("dummy","",colnames(dummy))
          mm=diag(ncol(dummy)); rownames(mm) <- colnames(mm) <- colnames(dummy)#namess2[i]
        }
        init2[[i]] <- list(dummy,mm)
      }
    }
  }
  # make a dataframe with the vectors and matrices provided by the user
  
  nre <- length(init2)
  Z <- init2[[length(init2)]][[1]]
  if(nre > 1){ # there's a structure
    strlist <- lapply(init2[1:(nre-1)], function(x){x[[2]]})
    if(length(strlist) >1){
      vcs <- do.call(function(...){kronecker(...,make.dimnames = TRUE)},strlist)
    }else{
      vcs <- strlist[[1]]
    }
  }
  if(nre==1){
    allzs <- matrix(1,nrow=nrow(Z),ncol=1); colnames(allzs) <- "u"
    vcs <- matrix(1,1,1); colnames(vcs) <- rownames(vcs) <- "u"
  }else{
    zs <- lapply(init2[1:(nre-1)], function(x){x[[1]]})
    allzs <- do.call(cbind,zs)
  }
  
  ## start creating the Z and K list
  Zup <- list()
  Kup <- list()
  typevc <- numeric()
  re_name <- character()
  counter <- 1
  for(i in 1:ncol(vcs)){ ## for each row
    for(j in 1:i){ ## for each column
      # print(paste(i,j))
      if(vcs[i,j] > 0){ ## to be estimated
        
        if(i==j){## var
          # commonlevs <- intersect(colnames(allzs),namz)
          # if(length(commonlevs) == 0){stop(paste("You may not be using a special variance structure in",paste(namess2,collapse = ","),"combination"),call. = FALSE)}
          namz <- strsplit(rownames(vcs)[i],":")[[1]]
          # print(matrix(apply(allzs[,namz],1,prod)))
          zz <- as.matrix(apply(as.matrix(allzs[,namz]),1,prod) * Z)
          if(is.null(Gu)){
            Gux <- diag(ncol(Z))
          }else{
            colnames(zz) <- gsub(ref_name,"",colnames(zz))
            checkg <- setdiff(colnames(zz),colnames(Gu))
            if(length(checkg)>0){
              stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
            }
            nameszz <- colnames(zz)
            Gux <- Gu[nameszz,nameszz]
          }
          Zup[[counter]] <- zz
          Kup[[counter]] <- Gux
          typevc[counter] <- 1
          re_name[counter] <- paste(rownames(vcs)[i],ref_name,sep=":")
          counter <- counter + 1
        }else{## cov
          namz1 <- strsplit(rownames(vcs)[i],":")[[1]]
          namz2 <- strsplit(colnames(vcs)[j],":")[[1]]
          z1 <- as.matrix(apply(as.matrix(allzs[,namz1]),1,prod) * Z)
          z2 <- as.matrix(apply(as.matrix(allzs[,namz2]),1,prod) * Z)
          if(is.null(Gu)){
            Gux <- diag(ncol(Z))
            Gu0 <- Gux*0
            Gu1 <- rbind(cbind(Gu0,Gux),cbind(Gux,Gu0))
          }else{
            colnames(z1) <- gsub(ref_name,"",colnames(z1))
            colnames(z2) <- gsub(ref_name,"",colnames(z2))
            checkg <- setdiff(colnames(z1),colnames(Gu))
            if(length(checkg)>0){
              stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
            }
            nameszz <- colnames(z1)
            Gu <- Gu[nameszz,nameszz]
            Gu0 <- Gux*0
            Gu1 <- rbind(cbind(Gu0,Gux),cbind(Gux,Gu0))
          }
          
          zz <- cbind(z1,z2)
          if(is.residual){ ## if residual we need to make Z square because we provide Zunits as the R
            zz <- zz%*% Gu1 %*% t(zz)
            Gu1 <- diag(ncol(zz))
          }
          Zup[[counter]] <- zz
          Kup[[counter]] <- Gu1
          typevc[counter] <- 2
          re_name[counter] <- paste(rownames(vcs)[i],colnames(vcs)[j],ref_name,sep=":")
          counter <- counter + 1
        }
      }
    }
  }
  
  if(is.null(Gtc)){
    if(!is.null(Gt)){
      if(is.list(Gt)){
        Gtc <- lapply(Gt, function(x){
          nt <- ncol(x)
          mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
          return(mm)
        })
      }else{
        nt <- ncol(Gt)
        mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
        Gtc <- mm
      }
    }
  }
  
  if(is.null(Gt)){ # user didn't provide Gt
    # Gt[lower.tri(Gt)] <- 0
    if(!is.null(Gtc)){ # user did provide Gtc so we need to complete them
      
      if(is.list(Gtc)){ ## if user provided a list
        
        if(is.residual){ftu <- 5}else{ftu <- 1}
        Gt <- lapply(Gtc,function(x){
          nt <- ncol(x)
          bnmm <- matrix(0.1,nt,nt)+diag(.05,nt)
          com <- (x/x); com[which(is.nan(com),arr.ind = TRUE)] <- 0
          if((ncol(bnmm) == ncol(com)) & (nrow(bnmm) == nrow(com)) ){ # random
            mm <- (bnmm*ftu)*com
          }else{mm <- bnmm}#fixed
        })
        
      }else{ # user provided a matrix
        
        nt <- ncol(Gtc)
        if(is.residual){
          # mm <- ( matrix(1,nt,nt) * 0 + 1) * 0.04977728 + diag(0.02488864, nt,nt)
          bnmm <- matrix(0.1,nt,nt)+diag(.05,nt)
          # print(Gtc)
          com <- (Gtc/Gtc); com[which(is.nan(com),arr.ind = TRUE)] <- 0
          if((ncol(bnmm) == ncol(com)) & (nrow(bnmm) == nrow(com)) ){ # random
            mm <- (bnmm*5)*com
          }else{mm <- bnmm}#fixed
        }else{
          bnmm <- matrix(0.1,nt,nt)+diag(.05,nt)
          # print(Gtc)
          com <- (Gtc/Gtc); com[which(is.nan(com),arr.ind = TRUE)] <- 0
          if((ncol(bnmm) == ncol(com)) & (nrow(bnmm) == nrow(com)) ){ # random
            mm <- bnmm*com
          }else{mm <- bnmm}#fixed
          
          # mm <- (matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt)
        }
        Gt <- mm
        
      }
      
    }
  }
  # S3$Gt <- Gt
  # S3$Gtc <- Gtc
  # S3$vcs <- vcs
  # Gtc <- lapply(Gtc,function(x){x[lower.tri(x)] <- 0; return(x)})
  S3 <- list(Z=Zup,K=Kup,Gt=Gt,Gtc=Gtc,typevc=typevc,re_name=re_name,vcs=vcs)
  return(S3)
}

spl2D <-  function(x.coord,y.coord,at,at.levels, type="PSANOVA", nseg = c(10,10), pord = c(2,2), degree = c(3,3), nest.div = c(1,1) ) {
  
  if(!is.numeric(x.coord)){
    stop("x.coord argument in spl2D() needs to be numeric.", call. = FALSE)
  }
  if(!is.numeric(y.coord)){
    stop("y.coord argument in spl2D() needs to be numeric.", call. = FALSE)
  }
  
  interpret.covarrubias.formula <-
    function(formula) {
      env <- environment(formula) 
      if(inherits(formula, "character"))          
        formula <- as.formula(formula)
      tf <- terms.formula(formula, specials = c("SAP", "PSANOVA"))
      terms <- attr(tf, "term.labels")
      nt <- length(terms)
      if(nt != 1)
        stop("Error in the specification of the spatial effect: only a sigle bidimensional function is allowed")
      
      res <- eval(parse(text = terms[1]), envir = env)
      res
    }
  
  bbase <-
    function(X., XL., XR., NDX., BDEG.) {
      # Function for B-spline basis
      dx <- (XR. - XL.)/NDX.
      knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
      P <- outer(X., knots, tpower, BDEG.)
      n <- dim(P)[2]
      D <- diff(diag(n), diff = BDEG. + 1) / (gamma(BDEG. + 1) * dx ^ BDEG.)
      B <- (-1) ^ (BDEG. + 1) * P %*% t(D)
      res <- list(B = B, knots = knots)
      res 
    }
  
  tpower <-
    function(x, t, p) {
      # Function for truncated p-th power function
      return((x - t) ^ p * (x > t))
    }
  
  Rten2 <-
    function(X1,X2) {
      one.1 <- matrix(1,1,ncol(X1))
      one.2 <- matrix(1,1,ncol(X2))
      kronecker(X1,one.2)*kronecker(one.1,X2)
    }
  
  MM.basis <-
    function (x, xl, xr, ndx, bdeg, pord, decom = 1) {
      Bb = bbase(x,xl,xr,ndx,bdeg)
      knots <- Bb$knots
      B = Bb$B
      m = ncol(B)
      n = nrow(B)
      D = diff(diag(m), differences=pord)
      P.svd = svd(crossprod(D))
      U.Z = (P.svd$u)[,1:(m-pord)] # eigenvectors
      d = (P.svd$d)[1:(m-pord)]  # eigenvalues
      Z = B%*%U.Z
      U.X = NULL
      if(decom == 1) {
        U.X = ((P.svd$u)[,-(1:(m-pord))])
        X = B%*%U.X
      } else if (decom == 2){
        X = NULL
        for(i in 0:(pord-1)){
          X = cbind(X,x^i)
        }
      } else if(decom == 3) {
        U.X = NULL
        for(i in 0:(pord-1)){
          U.X = cbind(U.X,knots[-c((1:pord),(length(knots)- pord + 1):length(knots))]^i)
        }
        X = B%*%U.X
      } else if(decom == 4) { # Wood's 2013
        X = B%*%((P.svd$u)[,-(1:(m-pord))])
        id.v <- rep(1, nrow(X))
        D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
        Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
        X <- X%*%Xf
        U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
      }
      list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots, U.X = U.X, U.Z = U.Z)
    }
  
  ####################
  ### if we want to use at and at.levels
  if(!missing(at)){
    col1 <- deparse(substitute(at))
    dat <- data.frame(x.coord, y.coord, at); colnames(dat) <- c("x.coord","y.coord",col1)
    by <- col1
    if(!missing(at.levels)){
      by.levels=at.levels
    }else{by.levels=NULL}
  }else{
    by=NULL
    by.levels=NULL
    dat <- data.frame(x.coord,y.coord); colnames(dat) <- c("x.coord","y.coord")
  }
  #######################
  
  x.coord <- "x.coord"
  y.coord <- "y.coord"
  
  if(is.null(by)){
    dat$FIELDINST <- "FIELD1"
    by="FIELDINST"
    dat[,by] <- as.factor(dat[,by])
    data0 <- split(dat, dat[,by])
    if(!is.null(by.levels)){
      keep <- names(data0)[which(names(data0) %in% by.levels)]
      
      if(length(keep)==0){stop("The by.levels provided were not found in your dataset.",call. = FALSE)}
      
      data0 <- data0[[keep]]
      if(length(keep)==1){data0 <- list(data0); names(data0) <- keep}
    }
  }else{
    check <- which(colnames(dat)==by)
    if(length(check)==0){stop("by argument not found in the dat provided", call. = FALSE)}else{
      
      missby <- which(is.na(dat[,by]))
      if(length(missby)>0){stop("We will split using the by argument and you have missing values in this column.\nPlease correct.", call. = FALSE)}
      
      dat[,by] <- as.factor(dat[,by])
      data0 <- split(dat, dat[,by])
      
      if(!is.null(by.levels)){
        keep <- names(data0)[which(names(data0) %in% by.levels)]
        
        if(length(keep)==0){stop("The by.levels provided were not found in your dataset.",call. = FALSE)}
        
        data0 <- data0[[keep]]
        if(length(keep)==1){data0 <- list(data0); names(data0) <- keep}
        #print(str(data0))
      }
      
    }
  }
  
  nasx <- which(is.na(dat[,x.coord]))
  nasy <- which(is.na(dat[,y.coord]))
  if(length(nasx) > 0 | length(nasy) >0){
    stop("x.coord and y.coord columns cannot have NA's", call. = FALSE)
  }
  #res <- interpret.covarrubias.formula(formula)
  
  ####
  #### now apply the same to all environments
  multires <- lapply(data0, function(data){
    
    
    x1 <- data[ ,x.coord]
    x2 <- data[ ,y.coord]
    
    #type = type
    
    MM1 = MM.basis(x1, min(x1), max(x1), nseg[1], degree[1], pord[1], 4)
    MM2 = MM.basis(x2, min(x2), max(x2), nseg[2], degree[2], pord[2], 4)
    
    X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
    X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B
    
    c1 = ncol(B1); c2 = ncol(B2)
    
    # Nested bases
    if(nest.div[1] == 1) {
      MM1n <- MM1
      Z1n <- Z1
      c1n <- c1
      d1n <- d1	
    } else {
      MM1n = MM.basis(x1, min(x1), max(x1), nseg[1]/nest.div[1], degree[1], pord[1], 4)
      Z1n <- MM1n$Z
      d1n <- MM1n$d
      c1n <-  ncol(MM1n$B)  					
    }
    if(nest.div[2] == 1) {
      MM2n <- MM2
      Z2n <- Z2
      c2n <- c2
      d2n <- d2	
    } else {
      MM2n = MM.basis(x2, min(x2), max(x2), nseg[2]/nest.div[2], degree[2], pord[2], 4)
      Z2n <- MM2n$Z
      d2n <- MM2n$d
      c2n <-  ncol(MM2n$B)  					
    }
    
    x.fixed <- y.fixed <- ""
    for(i in 0:(pord[1]-1)){
      if(i == 1) 
        x.fixed <- c(x.fixed, x.coord)
      else if( i > 1)
        x.fixed <- c(x.fixed, paste(x.coord, "^", i, sep = ""))
    }
    for(i in 0:(pord[2]-1)){
      if(i == 1) 
        y.fixed <- c(y.fixed, y.coord)
      else if( i > 1)
        y.fixed <- c(y.fixed, paste(y.coord, "^", i, sep = ""))
    }
    xy.fixed <- NULL
    for(i in 1:length(y.fixed)) {
      xy.fixed <- c(xy.fixed, paste(y.fixed[i], x.fixed, sep= ""))
    }
    xy.fixed <- xy.fixed[xy.fixed != ""]
    names.fixed <- xy.fixed
    
    smooth.comp <- paste("f(", x.coord,",", y.coord,")", sep = "")
    
    if(type == "SAP") {
      names.random <- paste(smooth.comp, c(x.coord, y.coord), sep = "|")				
      X = Rten2(X2, X1)		
      # Delete the intercept
      X <- X[,-1,drop = FALSE]
      Z = cbind(Rten2(X2, Z1), Rten2(Z2, X1), Rten2(Z2n, Z1n))
      
      dim.random <- c((c1 -pord[1])*pord[2] , (c2 - pord[2])*pord[1], (c1n - pord[1])*(c2n - pord[2]))		
      dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
      names(dim$fixed) <- names.fixed
      names(dim$random) <- paste(smooth.comp, "Global")
      
      # Variance/Covariance components
      g1u <- rep(1, pord[2])%x%d1
      g2u <- d2%x%rep(1, pord[1])
      g1b <- rep(1, c2n - pord[2])%x%d1n
      g2b <- d2n%x%rep(1, c1n - pord[1])
      
      g <- list()	
      g[[1]] <- c(g1u, rep(0, dim.random[2]), g1b)
      g[[2]] <- c(rep(0, dim.random[1]), g2u, g2b)
      
      names(g) <- names.random
      
    } else {		
      one1. <- X1[,1, drop = FALSE]
      one2. <- X2[,1, drop = FALSE]
      
      x1. <- X1[,-1, drop = FALSE]
      x2. <- X2[,-1, drop = FALSE]
      
      # Fixed and random matrices
      X <- Rten2(X2, X1)
      # Delete the intercept
      X <- X[,-1,drop = FALSE]
      Z <- cbind(Rten2(one2., Z1), Rten2(Z2, one1.), Rten2(x2., Z1), Rten2(Z2, x1.), Rten2(Z2n, Z1n))
      
      dim.random <- c((c1-pord[1]), (c2-pord[2]), (c1-pord[1])*(pord[2]-1), (c2-pord[2])*(pord[1]-1), (c1n-pord[2])*(c2n-pord[2]))
      
      # Variance/Covariance components		
      g1u <- d1
      g2u <- d2
      
      g1v <- rep(1, pord[2] - 1)%x%d1
      g2v <- d2%x%rep(1,pord[1] - 1)
      
      g1b <- rep(1, c2n - pord[2])%x%d1n
      g2b <- d2n%x%rep(1, c1n - pord[1])
      
      g <- list()
      
      if(type == "SAP.ANOVA") {
        g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
        g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
        g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, dim.random[4]), g1b)
        g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, g2b)
        
        names.random <- c(paste("f(", x.coord,")", sep = ""), paste("f(", y.coord,")", sep = ""), paste(smooth.comp, c(x.coord, y.coord), sep = "|"))			
        dim <- list(fixed = rep(1, ncol(X)), random = c(dim.random[1:2], sum(dim.random[-(1:2)])))		
        names(dim$fixed) <- names.fixed
        names(dim$random) <- c(names.random[1:2], paste(smooth.comp, "Global"))
        names(g) <- names.random
      } else {
        g[[1]] <- c(g1u, rep(0, sum(dim.random[2:5])))
        g[[2]] <- c(rep(0, dim.random[1]), g2u, rep(0, sum(dim.random[3:5])))
        g[[3]] <- c(rep(0, sum(dim.random[1:2])), g1v, rep(0, sum(dim.random[4:5])))
        g[[4]] <- c(rep(0, sum(dim.random[1:3])), g2v, rep(0, dim.random[5]))
        g[[5]] <- c(rep(0, sum(dim.random[1:4])), g1b + g2b)
        
        names.random <- c(paste("f(", x.coord,")", sep = ""), paste("f(", y.coord,")", sep = ""),
                          paste("f(", x.coord,"):", y.coord, sep = ""),
                          paste(x.coord,":f(", y.coord,")", sep = ""),
                          paste("f(", x.coord,"):f(", y.coord,")", sep = ""))
        
        dim <- list(fixed = rep(1, ncol(X)), random = dim.random)		
        names(dim$fixed) <- names.fixed
        names(dim$random) <- names.random
        names(g) <- names.random
      }		
    }
    colnames(X) <- names.fixed
    colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
    
    attr(dim$fixed, "random") <- attr(dim$fixed, "sparse") <- rep(FALSE, length(dim$fixed))
    attr(dim$fixed, "spatial") <- rep(TRUE, length(dim$fixed))
    
    attr(dim$random, "random") <- attr(dim$random, "spatial") <- rep(TRUE, length(dim$random)) 
    attr(dim$random, "sparse") <- rep(FALSE, length(dim$random))
    
    terms <- list()
    terms$MM <- list(MM1 = MM1, MM2 = MM2)
    terms$MMn <- list(MM1 = MM1n, MM2 = MM2n)
    #terms$terms.formula <- res
    
    # attr(terms, "term") <- smooth.comp
    
    # Initialize variance components
    init.var <- rep(1, length(g))
    
    res <- list(X = X, Z = Z, dim = dim, g = g, init.var = init.var)	
    M <- cbind(res$X,res$Z)
    
    return(M)
  })
  
  nrowss <- (unlist(lapply(multires,nrow)))
  nranges <- (unlist(lapply(multires,ncol)))
  
  names(multires) <- gsub(" ",".",names(multires))
  names(multires) <- gsub("#",".",names(multires))
  names(multires) <- gsub("-",".",names(multires))
  names(multires) <- gsub("/",".",names(multires))
  names(multires) <- gsub("%",".",names(multires))
  names(multires) <- gsub("\\(",".",names(multires))
  names(multires) <- gsub(")",".",names(multires))
  
  
  st <- 1 # for number of rows
  st2 <- 1 # for number of column
  end2 <- numeric() # for end of number of column
  dataflist <- list()
  #glist <- list()
  for(u in 1:length(multires)){
    prov <- multires[[u]]
    mu <- as.data.frame(matrix(0,nrow = sum(nrowss), ncol = ncol(prov)))
    colnames(mu) <- paste(names(multires)[u],colnames(prov), sep="_")
    
    nam <- paste("at",names(multires)[u],"2Dspl", sep="_")
    end <- as.numeric(unlist(st+(nrowss[u]-1)))
    
    mu[st:end,] <- prov
    dataflist[[nam]] <- as(as.matrix(mu), Class="sparseMatrix")
    st <- end+1
    ## for keeping track of the inits
    # end2 <- as.numeric(unlist(st2+(ncol(prov)-1)))
    # glist[[nam]] <- st2:end2
    # st2 <- end2+1
  }
  
  ## now build the last dataframe and adjust the glist
  # newdatspl <- as.data.frame(do.call(cbind,dataflist))
  # nn <- ncol(dat) # to add to the glist
  # glist <- lapply(glist, function(x){x+nn})
  # newdat <- data.frame(dat,newdatspl)
  
  ## now make the formula
  
  #funny <- paste(paste("grp(",names(dataflist),")",sep=""), collapse=" + ")
  
  ## important
  # newdat: is the neew data frame with original data and splines per location matrices
  # glist: is the argument to provide in group in asreml to indicate where each grouping starts and ends
  # funny: formula to add to your random formula
  dataflist <- lapply(dataflist,as.matrix)
  fin <- Reduce("+",dataflist)
  # fin <-dataflist#list(newdat=dataflist, funny=funny) # 
  return(fin)
}

bivariateRun <- function(model, n.core=1){
  
  args <- model[[length(model)]]
  
  if(!args$return.param){
    stop("The model provided needs to have the return.param argument set to TRUE. \nPlease read the documentation of the bivariateRun function carefully.\n", call. = FALSE)
  }
  
  response <- strsplit(as.character(args$fixed[2]), split = "[+]")[[1]]
  expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
  traits <- trimws(strsplit(expi(response),",")[[1]])
  
  combos <- expand.grid(traits,traits)
  combos <- combos[which(combos[,1] != combos[,2]),]; 
  combos <- combos[!duplicated(t(apply(combos, 1, sort))),];rownames(combos) <- NULL
  
  RHS <- as.character(args$fixed[3])
  it <- as.list(1:nrow(combos))
  
  cat(paste(nrow(combos), "bivariate models to be run\n"))
  
  model.results <- parallel::mclapply(it, 
                                      function(x) {
                                        # score.calc(M[ix.pheno, markers])
                                        ff <- as.formula(paste("cbind(",paste(as.vector(unlist(combos[x,])),collapse = ","),") ~", RHS))
                                        # do the fixed call
                                        args2 <- args; args2$fixed <- ff; args2$return.param <- FALSE
                                        # modify the random call for 2 traits
                                        p1 <- gsub("unsm\\([[:digit:]])","unsm(2)",as.character(args2$random))
                                        p1 <- gsub("diag\\([[:digit:]])","diag(2)",p1)
                                        p1 <- gsub("uncm\\([[:digit:]])","uncm(2)",p1)
                                        args2$random <- as.formula(paste(p1[1],paste(p1[-1],collapse = "+")))
                                        # modify the rcov call for 2 traits
                                        p1 <- gsub("unsm\\([[:digit:]])","unsm(2)",as.character(args2$rcov))
                                        p1 <- gsub("diag\\([[:digit:]])","diag(2)",p1)
                                        p1 <- gsub("uncm\\([[:digit:]])","uncm(2)",p1)
                                        args2$rcov <- as.formula(paste(p1[1],paste(p1[-1],collapse = "+")))
                                        
                                        gsub("[1-9]","k",as.character(args2$random))
                                        res0 <- do.call(mmer, args=args2)
                                        return(res0)
                                      }, 
                                      mc.cores = n.core)
  
  sigmas <- lapply(model.results, function(x){x$sigma})
  sigmas_scaled <- lapply(model.results, function(x){x$sigma_scaled})
  nre <- length(sigmas[[1]])
  namesre <- names(model.results[[1]]$sigma)
  sigmaslist <- list()
  sigmas_scaledlist <- list()
  for(i in 1:nre){
    mt <- matrix(0,length(traits),length(traits))
    rownames(mt) <- colnames(mt) <- traits
    mts <- mt
    sigmaprov <- lapply(sigmas, function(x){x[[i]]})
    sigmascaledprov <- lapply(sigmas_scaled, function(x){x[[i]]})
    for(j in 1:length(sigmaprov)){
      mt[colnames(sigmaprov[[j]]),colnames(sigmaprov[[j]])] <- sigmaprov[[j]]
      mts[colnames(sigmascaledprov[[j]]),colnames(sigmascaledprov[[j]])] <- sigmascaledprov[[j]]
    }
    sigmaslist[[namesre[i]]] <- mt
    sigmas_scaledlist[[namesre[i]]] <- mts
  }
  
  names(model.results) <- apply(combos,1,function(x){paste(x,collapse = "-")})
  # corlist <- lapply(sigmaslist,cov2cor)
  final <- list(sigmas=sigmaslist, sigmas_scaled=sigmas_scaledlist, models=model.results)
  return(final)
}

transformConstraints <- function(list0,value=1){
  ll <- lapply(list0, function(x){
    x[which(x != 0,arr.ind = TRUE)] <- x[which(x != 0,arr.ind = TRUE)] / x[which(x != 0,arr.ind = TRUE)]
    x <- x*value
    return(x)
  })
  return(ll)
}

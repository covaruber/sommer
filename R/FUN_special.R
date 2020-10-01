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
  ## keep track of factor variables
  myTypes <- unlist(lapply(init,class))
  init0 <- init
  ##
  init <- lapply(init, as.character)
  names <- as.character(substitute(list(...)))[-1L]
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
  attr(S3,"variables") <- names
  return(S3)
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

vs <- function(..., Gu=NULL, Gti=NULL, Gtc=NULL){
  
  ## ... list of structures to define the random effect
  ## Gu the known covariance matrix of the vs
  ## Gti the multitrait structure and constraints for it
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
  # certain random effects coming from spl2D(), leg(), and others may need some help to find the terms
  specialVariables <- unlist(lapply(init,function(x){(attributes(x)$variables)}))
  # print(namess2)
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
        # print(str(init[[i]]))
        # print(attributes(init[[i]])$variables)
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
            # colnames(zz) <- gsub(ref_name,"",colnames(zz)) ## why I wrote this?
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
            # colnames(z1) <- gsub(ref_name,"",colnames(z1)) ## why I wrote this?
            # colnames(z2) <- gsub(ref_name,"",colnames(z2)) ## why I wrote this?
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
    if(!is.null(Gti)){
      if(is.list(Gti)){
        Gtc <- lapply(Gti, function(x){
          nt <- ncol(x)
          mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
          return(mm)
        })
      }else{
        nt <- ncol(Gti)
        mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
        Gtc <- mm
      }
    }
  }
  
  if(is.null(Gti)){ # user didn't provide Gti
    # Gti[lower.tri(Gti)] <- 0
    if(!is.null(Gtc)){ # user did provide Gtc so we need to complete them
      
      if(is.list(Gtc)){ ## if user provided a list
        
        if(is.residual){ftu <- 5}else{ftu <- 1}
        Gti <- lapply(Gtc,function(x){
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
        Gti <- mm
        
      }
      
    }
  }
  # S3$Gti <- Gti
  # S3$Gtc <- Gtc
  # S3$vcs <- vcs
  # Gtc <- lapply(Gtc,function(x){x[lower.tri(x)] <- 0; return(x)})
  if(!is.null(specialVariables)){
    namess2 <- specialVariables
  }
  S3 <- list(Z=Zup,K=Kup,Gti=Gti,Gtc=Gtc,typevc=typevc,re_name=re_name,vcs=vcs, terms=namess2)
  return(S3)
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

# unsBLUP <- function(blups){
#   l <- unlist(lapply(blups,function(x){length(x[[1]])}))
#   lmin <- min(l); lmax <- max(l)
#   indexCov1 <- 1:lmin
#   indexCov2 <- (lmin+1):lmax
#   ntraits <- length(blups[[1]])
#   # blups follow the order of a lower triangula matrix
#   # (n*(n-1))/2 = l
#   # l*2 = n2 - n
#   # n2 = l*2 - n
#   n <- 1:100
#   possibilities <- ((n*(n-1))/2) + n
#   ntrue <- n[which(possibilities == length(l))]
#   ## index to know how to add them up
#   base <- matrix(NA,ntrue,ntrue)
#   base[lower.tri(base, diag=TRUE)] <- 1:length(l)
#   index <- which(!is.na(base), arr.ind = TRUE)
#   index <- index[order(index[,1]), ]
#   
#   
#   for(i in 1:ntrue){ # for each main blup
#     main <- which(index[,1] == i & index[,2] == i, arr.ind = TRUE)
#     cov1 <- which(index[,1] == i & index[,2] != i, arr.ind = TRUE)
#     cov2 <- which(index[,1] != i & index[,2] == i, arr.ind = TRUE)
#     for(itrait in 1:ntraits){
#       start <- blups[[main]][[itrait]]  
#       for(icov1 in cov1){
#         start <- start + blups[[icov1]][[itrait]][indexCov1]
#       }
#       for(icov2 in cov2){
#         start <- start + blups[[icov2]][[itrait]][indexCov2]
#       }
#       # store adjusted blup adding covariance effects in the same structure
#       blups[[main]][[itrait]] <- start
#     }
#   }
#   return(blups)
# }

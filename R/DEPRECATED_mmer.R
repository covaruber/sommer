# Legacy mmer API: retained unchanged pending deprecation.
# The historical covm below is superseded by CovarianceFactor-v2 covm in FUN_vsm.R.

# ---- FUN_special.R ----
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





# ---- FUN_methods_mmer.R ----
#### =========== ####
## SUMMARY FUNCTION mmer #
#### =========== ####
"summary.mmer" <- function(object, ...) {
  # message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
  replace.values <- function(Values,Search,Replace){
    dd0 <- data.frame(Values)
    vv <- which(Values%in%Search)
    dd <- data.frame(Search,Replace)
    rownames(dd) <- Search
    dd0[vv,"Values"] <- as.character(dd[Values[vv],"Replace"])
    return(dd0[,1])
  }
  
  if(!object$reshapeOutput){stop("summary function only works for reshaped output.", call. = FALSE)}
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

### >>>>>>>>
### vsr compiler function

vsr <- function(..., Gu=NULL, Gti=NULL, Gtc=NULL, reorderGu=TRUE, buildGu=TRUE){
  
  # message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  ## ... list of structures to define the random effect , e.g. init <- list(ds(M$data$FIELD),TP)
  ## Gu the known covariance matrix of the vs
  ## Gti the multitrait structure and constraints for it
  ## Gtc the initial values for the var-cov components
  
  init <- list(...)
  namess <- as.character(substitute(list(...)))[-1L]
  # expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  # expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=TRUE)}
  
  namess2 <- apply(data.frame(namess),1,function(x){
    return(paste(all.vars(as.formula(paste0("~",x))), collapse = ":"))
  })
  namess2 <- unique(unlist(namess2))
  # namess2 <- apply(data.frame(namess),1,function(x){
  #   newx <- expi2(x); if(length(newx)==0){newx<-""} 
  #   # an issue using predict() for a model of the form vs(us(leg(x)),y) led us to use expi2 in 09/14/2021. I am still not sure if there's consequences 
  #   newx <- gsub(",.*","",newx)
  #   return(newx)
  # })
  namess2[which(namess2 == "")] <- namess[which(namess2 == "")]
  # print(namess2)
  ref_name <- namess2[length(namess2)]
  # certain random effects coming from spl2D(), leg(), and others may need some help to find the terms
  specialVariables <- unlist(lapply(init,function(x){(attributes(x)$variables)}))
  # print(namess2)
  if("units" %in% namess2){
    is.residual =TRUE
  }else{is.residual=FALSE}
  ### get the data
  init2 <- list() # store the incidence matrices and relationship matrices for each random effect provided in the element "..."
  for(i in 1:length(init)){
    if(is.list(init[[i]])){ ## if it comes from a ds, us, cs function
      
      init2[[i]] <- init[[i]]
      
    }else{ # is a single vector with numbers or characters, ...
      
      if(is.matrix(init[[i]])){ # a mtrix is provided already so no need to create it
        if(buildGu){ # if user want Gu built
          mm=diag(ncol(init[[i]])); rownames(mm) <- colnames(mm) <- colnames(init[[i]])
        }else{ # if user don't want to build Gu most likely because is an rrBLUP model (millions of SNPs)
          mm=diag(1); rownames(mm) <- colnames(mm) <- colnames("dummy")
        }
        init2[[i]] <- list(x=init[[i]],mm)
      }else{ # is a vector
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
  
  nre <- length(init2) # number of covariates involved
  Z <- init2[[length(init2)]][[1]] # get incidence matrix from the main effect involved
  if(nre > 1){ # there's a covariance structure
    # we obtain vcs
    strlist <- lapply(init2[1:(nre-1)], function(x){x[[2]]})
    if(length(strlist) >1){ # there's a complicated structure e.g., vs(ds(a),ds(b),...,x)
      vcs <- do.call(function(...){kronecker(...,make.dimnames = TRUE)},strlist)
    }else{ # there's a simpler structure e.g., vs(ds(a),x)
      vcs <- strlist[[1]]
    }
    # we obtain z's of the covariates building the covariance structure e.g. in vs(ds(a),x) zs = model.matrix(~a)
    zs <- lapply(init2[1:(nre-1)], function(x){x[[1]]})
    allzs <- do.call(cbind,zs) # bind all the zs since we have multiple covariates building the covariance structure
  }else{
    # we obtain vcs
    vcs <- matrix(1,1,1); colnames(vcs) <- rownames(vcs) <- "u"
    # we obtain z's of the covariates building the covariance structure
    allzs <- matrix(1,nrow=nrow(Z),ncol=1); colnames(allzs) <- "u"
  }
  ###################################
  ## start creating the Z and K list
  Zup <- list() # store incidence matrices
  Kup <- list() # store relationship matrices between levels in Z
  typevc <- numeric() # store wheter is a variance (1) or covariance (2;allowed to be negative) component
  re_name <- character() # store the name of the random effect
  counter <- 1
  # print(vcs)
  for(i in 1:ncol(vcs)){ ## for each row
    for(j in 1:i){ ## for each column
      # print(paste(i,j))
      if(vcs[i,j] > 0){ ## to be estimated so build the incidence matrix
        
        if(i==j){## variance component (diagonal term in vcs) and either positive or fixed
          namz <- strsplit(rownames(vcs)[i],":")[[1]] # level name to assign together with main name
          # print(dim(allzs));print(dim(Z))
          zz <- as.matrix(apply(as.matrix(allzs[,namz]),1,prod) * Z) # take the column "namz" for the level and multiply by the provided main Z matrix 
          # image(as(zz, Class="sparseMatrix"))
          if(is.null(Gu)){
            
            if(!is.null(Gtc)){ # warning for possible mistake
              # print(Gtc)
              # print("a")
              # print(Gtc[i,j])
              if(vcs[i,j] == 2){ # if the user provides a covariance component as a variance component warn him he might not be providing the proper relationship matrix
                warning("You have provided an unconstrained variance component for a term with diagonal structure. \nA customized relationship structure may be needed.",call. = FALSE)
              }
            }
            
            if(buildGu){ # if user want Gu built
              Gux <- diag(ncol(Z))
            }else{ # if user don't want to build Gu most likely because is an rrBLUP model (millions of SNPs)
              Gux <- diag(1)
            }
            
          }else{ # if Gu is provided
            
            checkg <- setdiff(colnames(zz),colnames(Gu))
            if(length(checkg)>0){
              stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
            }
            checkg2 <- setdiff(colnames(Gu),colnames(zz))
            if(length(checkg2)>0){
              if(i==1){cat(paste0("Adding additional levels of Gu in the model matrix of '",ref_name,"' \n"))}
              added <- matrix(0, nrow = nrow(zz), ncol = length(checkg2)); colnames(added) <- checkg2
              zz <- cbind(zz,added)
            }
            nameszz <- colnames(zz)
            Gux <- Gu[nameszz,nameszz]
            if(!reorderGu){ # fix possible mistake
              # if the user provides a covariance component as a variance component do not rearrange the relationship matrix
              Gux <- Gu
            }
            
          }
          Zup[[counter]] <- zz
          Kup[[counter]] <- Gux
          typevc[counter] <- 1
          re_name[counter] <- paste(rownames(vcs)[i],ref_name,sep=":")
          counter <- counter + 1
        }else{## covariance component
          # print("cov")
          namz1 <- strsplit(rownames(vcs)[i],":")[[1]] # name of term1
          namz2 <- strsplit(colnames(vcs)[j],":")[[1]] # name of term2
          z1 <- as.matrix(apply(as.matrix(allzs[,namz1]),1,prod) * Z)
          z2 <- as.matrix(apply(as.matrix(allzs[,namz2]),1,prod) * Z)
          
          if(is.null(Gu)){
            
            if(buildGu){ # if user want Gu built
              Gux <- diag(ncol(Z))
            }else{ # if user don't want to build Gu most likely because is an rrBLUP model (millions of SNPs)
              Gux <- diag(1)
            }
            Gu0 <- Gux*0
            Gu1 <- rbind(cbind(Gu0,Gux),cbind(Gux,Gu0)) # image(as(Gu1, Class="sparseMatrix"))
          }else{
            
            checkg <- setdiff(colnames(z1),colnames(Gu))
            if(length(checkg)>0){
              stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
            }
            checkg2 <- setdiff(colnames(Gu),colnames(z1))
            if(length(checkg2)>0){
              if(i==1){cat(paste0("Adding additional levels of Gu in the model matrix of '",ref_name,"' \n"))}
              added <- matrix(0, nrow = nrow(z1), ncol = length(checkg2)); colnames(added) <- checkg2
              z1 <- cbind(z1,added)
            }
            
            checkg <- setdiff(colnames(z2),colnames(Gu))
            if(length(checkg)>0){
              stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
            }
            checkg2 <- setdiff(colnames(Gu),colnames(z2))
            if(length(checkg2)>0){
              if(i==1){cat(paste0("Adding additional levels of Gu in the model matrix of '",ref_name,"' \n"))}
              added <- matrix(0, nrow = nrow(z2), ncol = length(checkg2)); colnames(added) <- checkg2
              z2 <- cbind(z2,added)
            }
            
            nameszz <- colnames(z1)
            Gu <- Gu[nameszz,nameszz]
            Gu0 <- Gux*0
            Gu1 <- rbind(cbind(Gu0,Gux),cbind(Gux,Gu0))
          }
          
          zz <- cbind(z1,z2)
          if(is.residual){ ## if residual we need to make Z square because we provide Zunits as the R
            zz <- zz%*% Gu1 %*% t(zz)
            
            if(buildGu){ # if user want Gu built
              Gu1 <- diag(ncol(zz))
            }else{ # if user don't want to build Gu most likely because is an rrBLUP model (millions of SNPs)
              Gu1 <- diag(1)
            }
            
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
  # print(is.null(specialVariables));print(namess2)
  if(!is.null(specialVariables)){
    namess2 <- unique(c(namess2,specialVariables))
  }
  S3 <- list(Z=Zup,K=Kup,Gti=Gti,Gtc=Gtc,typevc=typevc,re_name=re_name,vcs=vcs, terms=namess2)
  return(S3)
}

vs <- function(..., Gu=NULL, Gti=NULL, Gtc=NULL, reorderGu=TRUE, buildGu=TRUE){
  #
  # message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  ## ... list of structures to define the random effect , e.g. init <- list(ds(M$data$FIELD),TP)
  ## Gu the known covariance matrix of the vs
  ## Gti the multitrait structure and constraints for it
  ## Gtc the initial values for the var-cov components
  
  init <- list(...)
  namess <- as.character(substitute(list(...)))[-1L]
  # expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
  # expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=TRUE)}
  
  namess2 <- apply(data.frame(namess),1,function(x){
    return(paste(all.vars(as.formula(paste0("~",x))), collapse = ":"))
  })
  namess2 <- unique(unlist(namess2))
  # namess2 <- apply(data.frame(namess),1,function(x){
  #   newx <- expi2(x); if(length(newx)==0){newx<-""} 
  #   # an issue using predict() for a model of the form vs(us(leg(x)),y) led us to use expi2 in 09/14/2021. I am still not sure if there's consequences 
  #   newx <- gsub(",.*","",newx)
  #   return(newx)
  # })
  namess2[which(namess2 == "")] <- namess[which(namess2 == "")]
  # print(namess2)
  ref_name <- namess2[length(namess2)]
  # certain random effects coming from spl2D(), leg(), and others may need some help to find the terms
  specialVariables <- unlist(lapply(init,function(x){(attributes(x)$variables)}))
  # print(namess2)
  if("units" %in% namess2){
    is.residual =TRUE
  }else{is.residual=FALSE}
  ### get the data
  init2 <- list() # store the incidence matrices and relationship matrices for each random effect provided in the element "..."
  for(i in 1:length(init)){
    if(is.list(init[[i]])){ ## if it comes from a ds, us, cs function
      
      init2[[i]] <- init[[i]]
      
    }else{ # is a single vector with numbers or characters, ...
      
      if(is.matrix(init[[i]])){ # a mtrix is provided already so no need to create it
        if(buildGu){ # if user want Gu built
          mm=diag(ncol(init[[i]])); rownames(mm) <- colnames(mm) <- colnames(init[[i]])
        }else{ # if user don't want to build Gu most likely because is an rrBLUP model (millions of SNPs)
          mm=diag(1); rownames(mm) <- colnames(mm) <- colnames("dummy")
        }
        init2[[i]] <- list(x=init[[i]],mm)
      }else{ # is a vector
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
  
  nre <- length(init2) # number of covariates involved
  Z <- init2[[length(init2)]][[1]] # get incidence matrix from the main effect involved
  if(nre > 1){ # there's a covariance structure
    # we obtain vcs
    strlist <- lapply(init2[1:(nre-1)], function(x){x[[2]]})
    if(length(strlist) >1){ # there's a complicated structure e.g., vs(ds(a),ds(b),...,x)
      vcs <- do.call(function(...){kronecker(...,make.dimnames = TRUE)},strlist)
    }else{ # there's a simpler structure e.g., vs(ds(a),x)
      vcs <- strlist[[1]]
    }
    # we obtain z's of the covariates building the covariance structure e.g. in vs(ds(a),x) zs = model.matrix(~a)
    zs <- lapply(init2[1:(nre-1)], function(x){x[[1]]})
    allzs <- do.call(cbind,zs) # bind all the zs since we have multiple covariates building the covariance structure
  }else{
    # we obtain vcs
    vcs <- matrix(1,1,1); colnames(vcs) <- rownames(vcs) <- "u"
    # we obtain z's of the covariates building the covariance structure
    allzs <- matrix(1,nrow=nrow(Z),ncol=1); colnames(allzs) <- "u"
  }
  ###################################
  ## start creating the Z and K list
  Zup <- list() # store incidence matrices
  Kup <- list() # store relationship matrices between levels in Z
  typevc <- numeric() # store wheter is a variance (1) or covariance (2;allowed to be negative) component
  re_name <- character() # store the name of the random effect
  counter <- 1
  # print(vcs)
  for(i in 1:ncol(vcs)){ ## for each row
    for(j in 1:i){ ## for each column
      # print(paste(i,j))
      if(vcs[i,j] > 0){ ## to be estimated so build the incidence matrix
        
        if(i==j){## variance component (diagonal term in vcs) and either positive or fixed
          namz <- strsplit(rownames(vcs)[i],":")[[1]] # level name to assign together with main name
          # print(dim(allzs));print(dim(Z))
          zz <- as.matrix(apply(as.matrix(allzs[,namz]),1,prod) * Z) # take the column "namz" for the level and multiply by the provided main Z matrix 
          # image(as(zz, Class="sparseMatrix"))
          if(is.null(Gu)){
            
            if(!is.null(Gtc)){ # warning for possible mistake
              # print(Gtc)
              # print("a")
              # print(Gtc[i,j])
              if(vcs[i,j] == 2){ # if the user provides a covariance component as a variance component warn him he might not be providing the proper relationship matrix
                warning("You have provided an unconstrained variance component for a term with diagonal structure. \nA customized relationship structure may be needed.",call. = FALSE)
              }
            }
            
            if(buildGu){ # if user want Gu built
              Gux <- diag(ncol(Z))
            }else{ # if user don't want to build Gu most likely because is an rrBLUP model (millions of SNPs)
              Gux <- diag(1)
            }
            
          }else{ # if Gu is provided
            
            checkg <- setdiff(colnames(zz),colnames(Gu))
            if(length(checkg)>0){
              stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
            }
            checkg2 <- setdiff(colnames(Gu),colnames(zz))
            if(length(checkg2)>0){
              if(i==1){cat(paste0("Adding additional levels of Gu in the model matrix of '",ref_name,"' \n"))}
              added <- matrix(0, nrow = nrow(zz), ncol = length(checkg2)); colnames(added) <- checkg2
              zz <- cbind(zz,added)
            }
            nameszz <- colnames(zz)
            Gux <- Gu[nameszz,nameszz]
            if(!reorderGu){ # fix possible mistake
              # if the user provides a covariance component as a variance component do not rearrange the relationship matrix
              Gux <- Gu
            }
            
          }
          Zup[[counter]] <- zz
          Kup[[counter]] <- Gux
          typevc[counter] <- 1
          re_name[counter] <- paste(rownames(vcs)[i],ref_name,sep=":")
          counter <- counter + 1
        }else{## covariance component
          # print("cov")
          namz1 <- strsplit(rownames(vcs)[i],":")[[1]] # name of term1
          namz2 <- strsplit(colnames(vcs)[j],":")[[1]] # name of term2
          z1 <- as.matrix(apply(as.matrix(allzs[,namz1]),1,prod) * Z)
          z2 <- as.matrix(apply(as.matrix(allzs[,namz2]),1,prod) * Z)
          
          if(is.null(Gu)){
            
            if(buildGu){ # if user want Gu built
              Gux <- diag(ncol(Z))
            }else{ # if user don't want to build Gu most likely because is an rrBLUP model (millions of SNPs)
              Gux <- diag(1)
            }
            Gu0 <- Gux*0
            Gu1 <- rbind(cbind(Gu0,Gux),cbind(Gux,Gu0)) # image(as(Gu1, Class="sparseMatrix"))
          }else{
            
            checkg <- setdiff(colnames(z1),colnames(Gu))
            if(length(checkg)>0){
              stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
            }
            checkg2 <- setdiff(colnames(Gu),colnames(z1))
            if(length(checkg2)>0){
              if(i==1){cat(paste0("Adding additional levels of Gu in the model matrix of '",ref_name,"' \n"))}
              added <- matrix(0, nrow = nrow(z1), ncol = length(checkg2)); colnames(added) <- checkg2
              z1 <- cbind(z1,added)
            }
            
            checkg <- setdiff(colnames(z2),colnames(Gu))
            if(length(checkg)>0){
              stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
            }
            checkg2 <- setdiff(colnames(Gu),colnames(z2))
            if(length(checkg2)>0){
              if(i==1){cat(paste0("Adding additional levels of Gu in the model matrix of '",ref_name,"' \n"))}
              added <- matrix(0, nrow = nrow(z2), ncol = length(checkg2)); colnames(added) <- checkg2
              z2 <- cbind(z2,added)
            }
            
            nameszz <- colnames(z1)
            Gu <- Gu[nameszz,nameszz]
            Gu0 <- Gux*0
            Gu1 <- rbind(cbind(Gu0,Gux),cbind(Gux,Gu0))
          }
          
          zz <- cbind(z1,z2)
          if(is.residual){ ## if residual we need to make Z square because we provide Zunits as the R
            zz <- zz%*% Gu1 %*% t(zz)
            
            if(buildGu){ # if user want Gu built
              Gu1 <- diag(ncol(zz))
            }else{ # if user don't want to build Gu most likely because is an rrBLUP model (millions of SNPs)
              Gu1 <- diag(1)
            }
            
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
  # print(is.null(specialVariables));print(namess2)
  if(!is.null(specialVariables)){
    namess2 <- unique(c(namess2,specialVariables))
  }
  S3 <- list(Z=Zup,K=Kup,Gti=Gti,Gtc=Gtc,typevc=typevc,re_name=re_name,vcs=vcs, terms=namess2)
  return(S3)
}

vcsExtract <- function(object){
  nre <- length(object$sigma)
  namesre <- names(object$sigma)
  vcs <- list()
  for(i in 1:nre){
    toextract <- which(object$constraints[[i]] > 0,arr.ind = TRUE)
    vcs[[i]] <- object$sigma[[i]][toextract]
    names1 <- apply(toextract,1,function(x){paste(colnames(object$sigma[[i]])[x[1]],colnames(object$sigma[[i]])[x[2]], sep="-")})
    names2 <- paste(namesre[i],names1, sep=".")
    names(vcs[[i]]) <- names2
  }
  vcs <- unlist(vcs)
  return(vcs)
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

##############
## special structures for mmer > vsr
atr <- function(x, levs){
  # message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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
  return(list(Z=dummy,thetaC=mm))
}
csr <- function(x,mm){
  # message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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
  return(list(Z=dummy,thetaC=mm))
}
dsr <- function(x){
  # message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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
  return(list(Z=dummy,thetaC=mm))
}
usr <- function(x){
  # message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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
  return(list(Z=dummy,thetaC=mm))
}


# ---- DEPRECATED_GWAS.R ----
GWAS <- function(fixed, random, rcov, data, weights, W,
                 nIters=20, tolParConvLL = 1e-03, tolParInv = 1e-06,
                 init=NULL, constraints=NULL, method="NR",
                 getPEV=TRUE,
                 naMethodX="exclude",
                 naMethodY="exclude",
                 returnParam=FALSE,
                 dateWarning=TRUE,date.warning=TRUE,
                 verbose=FALSE,
                 stepWeight=NULL, emWeight=NULL,
                 M=NULL, gTerm=NULL, n.PC = 0, min.MAF = 0.05,
                 P3D = TRUE){

  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
  if(is.null(gTerm)){
    stop("Please provide the name of the genetic term in the model in the gTerm argument", call. = FALSE)
  }

  if(length(which(is.na(M))) > 0){
    stop("Please provide an imputed marker matrix (M).", call. = FALSE)
  }
  make.full <- function(X) {
    svd.X <- svd(X)
    r <- max(which(svd.X$d > 1e-08))
    return(as.matrix(svd.X$u[, 1:r]))
  }

  dwToUse <- c(dateWarning,date.warning)
  v<-which(!dwToUse)
  if(length(v) == 0){v<-1}
  # print("hey")
  ## return all parameters for a mixed model (NOT ACTUAL FITTING)
  res <- mmer(fixed=fixed, random=random, rcov=rcov, data=data, W=W,  # silence weights and W when testing
              nIters=nIters, tolParConvLL=tolParConvLL, tolParInv=tolParInv,
              init=init, constraints=constraints, method=method,
              getPEV=getPEV,
              naMethodX=naMethodX,
              naMethodY=naMethodY,
              returnParam=TRUE,
              dateWarning=dwToUse[v],
              verbose=verbose,
              stepWeight=stepWeight, emWeight=emWeight
              )
  # print("hey")
  # print(str(res))

  if(returnParam){ # if user only wants the initial parameters and objects
    lastmodel <- res
  }else{  # if actual fitting is required
    ## fit initial models if P3D
    if(P3D){ # same VC for all markers
      lastmodel <- .Call("_sommer_MNR",PACKAGE = "sommer",res$yvar, res$X,res$Gx,
                         res$Z,res$K,res$R,res$GES,res$GESI,res$W, res$isInvW,
                         res$nIters,res$tolParConvLL,res$tolParInv, res$selected,res$getPEV,res$verbose,
                         TRUE, res$stepWeight, res$emWeight)
      Vinv <- lastmodel$Vi
      ## get names of random effects
      re_names <- unlist(res$re_names)
      re_names <- gsub("\\(Intercept):","",re_names)
      gTermi <- which(re_names %in% gTerm)
      if(length(gTermi) < 1){
        stop("No match between the Gterm and the random effects present.")
      }
      ## get input matrices
      Y <- scale(res$yvar)
      Z <- do.call(cbind,res$Z[gTermi])

      if (n.PC > 0) {
        Kb <- do.call(adiag1,res$K[gTermi])
        eig.vec <- eigen(Kb)$vectors
        X <- do.call(cbind,res$X)
        zbeig <- Z %*% eig.vec[,1:n.PC]
        X <- make.full(cbind(X,zbeig))
      }else {
        X <- do.call(cbind,res$X)
        X <- make.full(X)
      }

      if(nrow(M) != ncol(Z)){
        stop(paste("Marker matrix M needs to have same numbers of rows(",nrow(M),") than columns of the gTerm incidence matrix(",ncol(Z),")."),call. = FALSE)
      }

      if(length(which(rownames(M) != colnames(Z)))){
        M <- M[colnames(Z),]
      }
      m <- ncol(M); colnamesM <- colnames(M)
      ######################
      ## scorecalc function
      ######################
      cat(red("Performing GWAS evaluation\n"))
      results <- .Call("_sommer_gwasForLoop",PACKAGE = "sommer",
                         M,Y,as.matrix(Z),X,Vinv,min.MAF,TRUE
      )
      preScores <- matrix(results[, , 1], nrow(results), ncol(results))
      bs <- matrix(results[, , 2], nrow(results), ncol(results))
      bs.se <- matrix(results[, , 3], nrow(results), ncol(results))

      v2 <- length(Y) - ((ncol(X)+1)*ncol(Y)) # ncol(XZMi)
      pvals <- pbeta(preScores, v2/2, 1/2)
      scores <- -log10(pvals)
      ########################
      rownames(scores) <- colnames(M)
      colnames(scores) <- colnames(Y)
      
      lastmodel$effects <- bs
      lastmodel$effects.se <- bs.se
      lastmodel$pvals <- pvals
      lastmodel$scores <- scores
      lastmodel$shape1 <- v2/2
      lastmodel$shape2 <- 1/2
      lastmodel$method <- method
      lastmodel$constraints <- res[[8]]
      class(lastmodel)<-c("mmergwas")
    }else{ # if different variance for each marker

      ## get names of random effects
      re_names <- unlist(res$re_names)
      re_names <- gsub("\\(Intercept):","",re_names)
      gTermi <- which(re_names %in% gTerm)
      if(length(gTermi) < 1){
        stop("No match between the Gterm and the random effects present.")
      }
      Y <- scale(res$yvar)
      # print(head(Y))
      Z <- do.call(cbind,res$Z[gTermi])
      if(nrow(M) != ncol(Z)){
        stop(paste("Marker matrix M needs to have same numbers of rows(",nrow(M),") than columns of the gTerm incidence matrix(",ncol(Z),")."),call. = FALSE)
      }

      if (n.PC > 0) {
        Kb <- do.call(adiag1,res$K[gTermi])
        eig.vec <- eigen(Kb)$vectors
        X <- do.call(cbind,res$X)
        zbeig <- Z %*% eig.vec[,1:n.PC]
        X <- make.full(cbind(X,zbeig))
      }else {
        X <- do.call(cbind,res$X)
        X <- make.full(X)
      }
      #
      preScores <- bs <- bs.se <- list()
      for(iMarker in 1:ncol(M)){
        # print(iMarker)
        mi <- M[,iMarker]
        if(var(mi) > 0){ # if var > 0
          Xmi <- cbind(X,mi)
          lastmodel <- .Call("_sommer_MNR",PACKAGE = "sommer",res$yvar, list(Xmi),res$Gx,
                             res$Z,res$K,res$R,res$GES,res$GESI,res$W, res$isInvW,
                             res$nIters,res$tolParConvLL,res$tolParInv, res$selected,res$getPEV,FALSE,
                             TRUE, res$stepWeight, res$emWeight)
          b = lastmodel$Beta
          var.b = diag(as.matrix(lastmodel$VarBeta))
          se.b = sqrt(abs(var.b))
          t.val  = b/se.b
          # print(t.val)
          preScores[[iMarker]] = 1 - pnorm(t.val[(nrow(t.val)-ncol(Y)+1):nrow(t.val),])
          #print(b); print(se.b)
          bs[[iMarker]] = b[(nrow(b)-ncol(Y)+1):nrow(b),] # last fixed effect where we put the marker
          bs.se[[iMarker]] = se.b[(length(se.b)-ncol(Y)+1):length(se.b)] # matching standard error
        }else{ # if var == 0
          preScores[[iMarker]] <- rep(1,ncol(Y))
          bs[[iMarker]] = rep(0,ncol(Y)) # effect
          bs.se[[iMarker]] = rep(Inf,ncol(Y))
        }
      } # for loop for each marker
      preScores <- do.call(rbind,preScores)
      bs <- do.call(rbind, bs)
      bs.se <- do.call(rbind, bs.se)
      scores <- apply(preScores, 2, function(x){-log10(x)})
      
      # scores <- as.matrix(-log10(preScores))
      
      rownames(scores) <- colnames(M) # marker names
      # print(head(scores));print(colnames(Y))
      colnames(scores) <- colnames(Y) # trait names
     
      lastmodel$effects <- bs
      lastmodel$effects.se <- bs.se
      lastmodel$scores <- scores
      lastmodel$pvals <- preScores
      lastmodel$method <- method
      lastmodel$constraints <- res[[8]]
      class(lastmodel)<-c("mmergwas")

    }

  } # if actual fitting is required

  return(lastmodel)
}

# ---- mmer.R ----
mmer <- function(fixed, random, rcov, data, weights, W,
                 nIters=20, tolParConvLL = 1e-03, tolParInv = 1e-06,
                 init=NULL, constraints=NULL, method="NR",
                 getPEV=TRUE,
                 naMethodX="exclude",
                 naMethodY="exclude",
                 returnParam=FALSE,
                 dateWarning=TRUE, date.warning=TRUE,
                 verbose=TRUE,reshapeOutput=TRUE,
                 stepWeight=NULL, emWeight=NULL, 
                 contrasts=NULL){

  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions. This function will be no longer maintained.")
  dateWarning=FALSE
  desc <- utils::packageDescription("sommer")
  my.date <- as.Date(desc$Date)+90
  # my.date <- "2026-06-01"
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
  }else{nodata=FALSE}
  data <- as.data.frame(data)
  if(missing(rcov)){
    rcov = as.formula("~units")
  }

  #################
  ## do the needed for naMethodY and naMethodX
  dataor <- data
  provdat <- subdata(data, fixed=fixed, na.method.Y = naMethodY,na.method.X=naMethodX)
  # print(str(provdat))
  data <- provdat$datar
  nonMissing <- provdat$good
  # print(length(nonMissing))
  # randomization <- sample(1:nrow(data))
  # data <- data[randomization,]
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
  yvar <- as.matrix(model.response(mfna))
  # print(head(yvar))
  nt <- ncol(yvar)
  if(nt==1){colnames(yvar) <- response}
  # print(colnames(yvar))
  #################
  ## get Zs and Ks
  termsR <- list()
  termsRN <- numeric()
  if(!missing(random)){
    yuyu <- strsplit(as.character(random[2]), split = "[+]")[[1]]
    rtermss <- apply(data.frame(yuyu),1,function(x){
      strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
    })
    # print(rtermss)
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
    for(u in 1:length(rtermss)){ # for each random effect

      checkvs <- intersect(all.names(as.formula(paste0("~",rtermss[u]))),c("vs","vsr","gvsr","spl2Da","spl2Db")) #
      # print(checkvs)
      if(length(checkvs)>0){ ## if this term is a variance structure
        # print(rtermss[u])
        # print(data)
        ff <- eval(parse(text = rtermss[u]),data,parent.frame())# envir = data)
        # print(nrow(ff$Z[[1]]))
        # print(length(provdat$good))
        # print((ff$Z[[1]]))
        if(nrow(ff$Z[[1]]) != length(provdat$good)){
          ## if the incidence matrix is different than the size of good very likely the user provided a matrix
          ff$Z <- lapply(ff$Z,function(xxx){as.matrix(xxx[provdat$good,])})
        }
        # print(rtermss[u])
        # print(ff$Gtc)
        re_namel1[[u]] <- ff$re_name # store the levels
        termsR[[u]] <- ff$terms # store the variables used for this variance structure
        termsRN[u] <- length(ff$Z)
        zs[[u]] <- ff$Z
        ks[[u]] <- ff$K
        if(is.null(ff$Gti)){ ## initial vc if user don't provide them
          mml <- list()
          for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
            if(ff$typevc[k] == 2){div=2}else{div=1}## divisor
            # mml[[k]] <- (( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt))/div
            mml[[k]] <- bnmm/div
          }
          ges[[u]] <- mml
        }else{ ## user provides initial values
          if(is.list(ff$Gti)){ ## user provides a list of initial values
            ges[[u]] <- ff$Gti
          }else{
            ges[[u]] <- rep(list(ff$Gti),length(ff$Z))
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
        re_namel1[[u]] <- list(rtermss[u]) # store the levels
        termsR[[u]] <- rtermss[u] # store the varibles used
        termsRN[u] <- 1
        ges[[u]] <- list(bnmm)
        mm <- matrix(1,nt,nt); mm[lower.tri(mm)] <- 0; mm[upper.tri(mm)] <- 2
        gesI[[u]] <- list(mm)
        # }
      }
    }
    Z <- unlist(zs,recursive=FALSE)
    K <- unlist(ks,recursive=FALSE)
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
  termsE <- list()
  termsEN <- numeric()
  bnmm <- matrix(0.1,nt,nt)+diag(.05,nt)
  for(u in 1:length(rcovtermss)){

    checkvs <- intersect(all.names(as.formula(paste0("~",rcovtermss[u]))),c("vs","vsr","gvsr"))

    if(length(checkvs)>0){ ## if this term is a variance structure
      ff <- eval(parse(text = rcovtermss[u]),data,parent.frame())# envir = data)
      rs[[u]] <- ff$Z
      re_namel2[[u]] <- ff$re_name
      termsE[[u]] <- ff$terms
      termsEN[u] <- length(ff$Z)
      if(is.null(ff$Gti)){ ## initial vc if user don't provide them
        mml <- list()
        for(k in 1:length(ff$Z)){ ## the diagonal of the matrix should be 1 is vc and 2 if cov
          if(ff$typevc[k] == 2){div=2}else{div=1}## divisor
          # mml[[k]] <- (( matrix(1,nt,nt) * 0 + 1) * 0.04977728 + diag(0.02488864, nt,nt))/div
          mml[[k]] <- (bnmm/div)*5
        }
        gesr[[u]] <-  mml
      }else{gesr[[u]] <- rep(list(ff$Gti),length(ff$Z))}
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
      termsE[[u]] <- rcovtermss[u]
      termsEN[u] <- 1
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

  # Expand the fixed terms to handle interactions
  fixedtermss <- attr(terms(fixed),"term.labels")

  # This will handle cases when user enters 1 anywhere in the fixed effects
  if(attr(terms(fixed),"intercept")==1){ # there should be intercept
    useinter <- TRUE
    fixedtermss <- c(1, fixedtermss)
    # }else if(test1 > 0 & test2 == 0){# there should be intercept as well
    # useinter <- TRUE
    # This will handle cases when user enters 1 anywhere in the fixed effects
  }else{ # there's no intercept
    useinter <- FALSE
    fixedtermss <- c(-1, fixedtermss)
  }
  # print(useinter)
  vsterms <- sort(unique(c(grep("vsr\\(",fixedtermss),grep("vs\\(",fixedtermss))), decreasing = FALSE)
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

  data$`1` <-1
  newfixed <- as.formula(paste("~",paste(fixedtermss,collapse = "+")))
  mf <- try(model.frame(newfixed, data = data, na.action = na.pass), silent = TRUE)
  mf <- eval(mf, parent.frame())
  baseX <- model.matrix(newfixed, mf, contrasts.arg=contrasts)
  # print(head(baseX))
  if(length(vsterms) > 0){baseX <- cbind(baseX,addxs)}
  if(as.double(nrow(baseX)) * as.double(ncol(baseX)) > 2147483647) {
    qr <- qr(baseX, LAPACK = TRUE)
  } else {
    qr <- qr(baseX)
  }
  keepx <- qr$pivot[1:qr$rank]
  nameskeepx <- colnames(baseX)[keepx]
  colsdropped <- ncol(baseX) - length(nameskeepx)
  baseX <- matrix(baseX[, keepx],nrow(yvar),qr$rank)
  colnames(baseX) <- nameskeepx
  if(colsdropped > 0){
    warning(paste("fixed-effect model matrix is rank deficient so dropping",colsdropped,"columns / coefficients\n"),call. = FALSE)
  }
  # print(colnames(baseX))
  xs <- list()
  gesf <- list()
  gesIf <- list()
  allfixedterms <- c(fixedtermss, fixedvsterms)
  for(u in 1:length(allfixedterms)){
    checkvs <- sort(unique(c(grep("vsr\\(",allfixedterms[u]),grep("vs\\(",allfixedterms[u]))), decreasing = FALSE)
    # checkvs <- grep("vsr\\(",allfixedterms[u])
    if(length(checkvs)>0){ ## if this term is a variance structure
      ffx <- eval(parse(text = allfixedterms[u]),data,parent.frame())# envir = data)
      findlevs <- colnames(ffx$Z[[1]])
      findlevs2 <- which(colnames(baseX) %in% findlevs)
      xpp <- as.matrix(baseX[, findlevs2])
      colnames(xpp) <- colnames(baseX)[findlevs2]
      # print(ff$Gtc)
      xs[[u]] <- list(xpp)
      if(is.null(ffx$Gti)){ ## initial vc if user don't provide them
        gesf[[u]] <- list(( matrix(1,nt,nt) * 0 + 1) * 0.1 + diag(0.05, nt))
      }else{gesf[[u]] <- rep(list(ffx$Gti),length(ffx$Z))}

      if(is.null(ffx$Gtc)){ ## contraints if user don't provide them
        mm <- diag(1,nt,nt); mm[lower.tri(mm)] <- 0; #mm[upper.tri(mm)] <- 2
        gesIf[[u]] <- list(mm)
      }else{gesIf[[u]] <- rep(list(ffx$Gtc),length(ffx$Z))}
    }else{ ## if is a normal term

      if(allfixedterms[u] == "1"){
        findlevs <- colnames(model.matrix(as.formula(paste("~",allfixedterms[u])), data=data, contrasts.arg=contrasts))
      }else{
        findlevsM <- colnames(model.matrix(as.formula(paste("~",allfixedterms[u],"-1")), data=data, contrasts.arg=contrasts))
        findlevsP <- setdiff( colnames(model.matrix(as.formula(paste("~",allfixedterms[u])), data=data, contrasts.arg=contrasts)) , "(Intercept)" )
        findlevs <- unique(c(findlevsM,findlevsP))
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
    # ws<- coco
    isInvW=TRUE
    WW <- diag(1/sqrt(coco)) # W is already squared and inverted
  }else{
    isInvW=TRUE
    ws <- rep(1,nrow(yvar)) # we build W as a diagonal
    WW <- diag(ws) #
  }
  if(!missing(W)){
    isInvW=FALSE
    WW <- W # user has provided W and needs to be squared and inverted
  }

  if(is.null(stepWeight)){
    stepWeight <- rep(0.9,nIters); stepWeight[1:2] <- c(0.5,0.7)
  }

  if(is.null(emWeight)){
    emWeight <- rep(0, nIters)
  }

  #################
  ## subset data
  if(method == "NR"){
    selected <- FALSE
  }else{selected <- TRUE}
  #################
  ## provide arguments to MNR
  if(!missing(random)){
    Z <- lapply(Z,function(x){as(as(as( x,  "dMatrix"), "generalMatrix"), "CsparseMatrix")})
  }else{Z <- list();K <- list();random=NULL}
  R <- lapply(R,function(x){as(as(as( x,  "dMatrix"), "generalMatrix"), "CsparseMatrix")})
  if(!is.null(init)){GES <- init}
  if(!is.null(constraints)){GESI <- constraints}
  re_names <- c(re_namel1,re_namel2)
  # print("good")
  if(returnParam){
    good <- provdat$good
    args <- list(fixed=fixed, random=random, rcov=rcov, data=data, weights=substitute(weights),
                 nIters=nIters, tolParConvLL=tolParConvLL, tolParInv=tolParInv,
                 init=init, constraints=constraints, method=method,
                 getPEV=getPEV,
                 naMethodX=naMethodX,
                 naMethodY=naMethodY,
                 returnParam=returnParam,
                 dateWarning=dateWarning,
                 verbose=verbose,reshapeOutput=reshapeOutput)

    res <- list(yvar=yvar, X=X,Gx=Gx,Z=Z,K=K,R=R,GES=GES,GESI=GESI, W=WW, isInvW=isInvW,
                nIters=nIters, tolParConvLL=tolParConvLL, tolParInv=tolParInv,
                selected=selected,getPEV=getPEV,verbose=verbose, retscaled=FALSE,
                re_names=re_names,good=good,fixedtermss=fixedtermss,args=args, stepWeight=stepWeight,
                emWeight=emWeight # emupdate=emupdate,
    )
  }else{

    res <- .Call("_sommer_MNR",PACKAGE = "sommer",yvar, X,Gx,Z,K,R,GES,GESI, WW, isInvW,
                 nIters, tolParConvLL, tolParInv,
                 selected,getPEV,verbose, FALSE, stepWeight, emWeight) # emupdate,

    # res <- MNR(yvar, X,Gx,Z,K,R,GES,GESI, W,isInvW,
    #              nIters, tolParConvLL, tolParInv,
    #              selected,getPEV,verbose, FALSE)

    if(length(res) > 0){
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
                       naMethodY=naMethodY,naMethodX=naMethodX)
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
      if(reshapeOutput){
        res <- reshape_mmer(res,namelist)
      }
      res$constraints <- GESI
      res$constraintsF <- Gx
      res$data <- data#dataor
      res$dataOriginal <- dataor
      res$terms <- list(response=list(colnames(yvar)),fixed=list(allfixedterms), random=termsR, rcov=termsE)
      res$termsN <- list(response=ncol(yvar),fixed=length(X), random=termsRN, rcov=termsEN)
      res$xEffectsN <- unlist(lapply(X,ncol))
      res$zEffectsN <- unlist(lapply(Z,ncol))
      fixedTerms <- unlist(res$terms$fixed)
      randomTerms <- unlist(res$terms$random)
      Dtable <- data.frame(type=c(rep("fixed",length(fixedTerms)), rep("random",length(randomTerms))),
                           term=c(fixedTerms, randomTerms),
                           include=FALSE, average=FALSE)
      res$Dtable <- Dtable
      res$W <- cbind(do.call(cbind,X),do.call(cbind,Z))
      # res$K <- K
      if(reshapeOutput){
        colnames(res$sigmaSE) <- rownames(res$sigmaSE) <- names(res$sigmaVector)
        res$sigmaVector <- vcsExtract(res)
      }
      res$reshapeOutput <- reshapeOutput

      class(res)<-c("mmer")
    }
  }
  # print("good")
  return(res)
}

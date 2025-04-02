#### =========== ####
## SUMMARY FUNCTION mmer #
#### =========== ####
"summary.mmer" <- function(object, ...) {
  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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

vsr <- function(..., Gu=NULL, Gti=NULL, Gtc=NULL, reorderGu=TRUE, buildGu=TRUE){
  
  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
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
  
  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
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

# myformula <- function(x){
#   expi <- function(j){gsub("[\\(\\)]", "", regmatches(j, gregexpr("\\(.*?\\)", j))[[1]])}
#   expi2 <- function(x){gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", x, perl=T)}
#   yuyuf <- strsplit(as.character(x[3]), split = "[+]")[[1]]
#   termss <- apply(data.frame(yuyuf),1,function(x){
#     strsplit(as.character((as.formula(paste("~",x)))[2]), split = "[+]")[[1]]
#   })
#   newtermss <- apply(data.frame(yuyuf),1,function(y){
#     newy <- expi(y)
#     if(length(newy) > 0){
#       newy <- gsub(",.*","",newy)
#     }else{newy <- y}
#     return(newy)
#   })
#   resp <- strsplit(as.character(x[2]), split = "[+]")[[1]]
#   newx <- paste(resp, "~",paste(newtermss,collapse = "+"))
#   return(newx)
# }

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
## VS structures
atr <- function(x, levs){
  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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
  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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
  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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
  message("This function has been deprecated. Please start using 'mmes' and its auxiliary functions (e.g., 'vsm', 'usm', 'dsm', 'ism', etc.). This function will be no longer maintained.")
  
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


vsc <- function(..., Gu=NULL, buildGu=TRUE, meN=1, meTheta=NULL, meThetaC=NULL, sp=FALSE){

  ## ... list of structures to define the random effect , e.g. init <- list(ds(M$data$FIELD),TP)
  ## Gu the known covariance matrix of the vs

  init <- list(...) #  e.g. init <- list(usc(data$Env),isc(data$Name)) | init <- list(dsc(data$YEAR),isc(data$units))

  namess <- as.character(substitute(list(...)))[-1L] # namess <- c("Env","Name") | namess <- c("YEAR","units")
  namess2 <- apply(data.frame(namess),1,function(x){
    return(all.vars(as.formula(paste0("~",x))))
  })

  ## let's test that user provided all terms encapsulated in a structure
  listLength <- length(init)
  if(listLength == 1){ # very simple structure
    myTest <- is.list(init[[1]])
    if(myTest){
      whichElemBad <- numeric()
    }else{whichElemBad <- 1}
  }else{ # there's more than one term
    lengthElem <- unlist(lapply(init,length))
    whichElemBad <- which(lengthElem < 3) ## which are not encapsulated
  }
  if(length(whichElemBad) > 0){
    badd <- paste(namess2[whichElemBad],collapse = ",")
    stop(paste0("Term(s): '",badd,"' in the vsc() function are not encapsulated in a structure function. Please correct [for example, using vsc(isc(",badd,")), vsc(dsc(",badd,")), vsc(usc(",badd,")), vsc(atc(",badd,")), vsc(csc(",badd,")), etc.]."),call. = FALSE)
  }

  ## extract names of variables and collpase as interaction

  namess2 <- as.vector(unique(unlist(namess2))) # remove repeats if exist
  namess2[which(namess2 == "")] <- namess[which(namess2 == "")] # remove empties if exist
  # extract the name of the main variable
  ref_name <- namess2[length(namess2)]
  # certain random effects coming from spl2D(), leg(), and others may need some help to find the terms
  specialVariables <- unlist(lapply(init,function(x){(attributes(x)$variables)}))
  #
  if("units" %in% namess2){ # if residual
    is.residual =TRUE
  }else{is.residual=FALSE}

  ### get the data
  if(length(init)>meN){ # if there is a covariance structure not only isc
    for(i in 1:(length(init)-meN)){ # for all terms prior to the main effect(s) keep nme in mind
      if(i==1){ # if is the first term
        theta <- init[[i]]$theta
        thetaC <- init[[i]]$thetaC
        Z0 <- init[[i]]$Z
      }else{ # after
        mm = init[[i]]$theta; mm=mm/mm; mm[which(is.nan(mm))]=0
        theta <- kronecker(theta,mm)
        thetaC <- kronecker(thetaC,init[[i]]$thetaC)
        provZ <- init[[i]]$Z
        provZlist <- list()
        for(j in 1:ncol(provZ)){
          # print(colnames(provZ))
          provZiCol <- provZ[,j] %*% Matrix(1,1,ncol(Z0))
          provZlist[[j]] <- Z0 * provZiCol
          colnames(provZlist[[j]]) <- paste(colnames(Z0),colnames(provZ)[j],sep=":")
        }
        Z0 <- cbind(Z0, do.call(cbind,provZlist))
      }
    }
    pasteNames=TRUE
  }else{ # if is only isc
    theta <- init[[1]]$theta
    thetaC <- init[[1]]$thetaC
    Z0 <- Matrix(1,nrow(init[[1]]$Z),1)
    pasteNames=FALSE
  }
  #######################################
  ## check covariance matrices
  `%!in%` = Negate(`%in%`)
  if(is.null(Gu)){
    x <- data.frame(d=as.factor(1:ncol(init[[length(init)]]$Z)))
    if(nrow(x) == 1){ # is a numeric variables
      Gu <- matrix(1);
      colnames(Gu) <- rownames(Gu) <- colnames(init[[length(init)]]$Z)
      Gu <- as(Gu, Class="dgCMatrix"); 
    }else{ # there's levels since it is a factor model
      Gu <- sparse.model.matrix(~d-1, x)
      colnames(Gu) <- rownames(Gu) <- colnames(init[[length(init)]]$Z)
    }
  }else{
    if (!inherits(Gu, "dgCMatrix")){
      stop("Gu matrix is not of class dgCMatrix. Please correct \n", call. = TRUE )
    }
  }
  #############################
  ######################################
  ## now build the Z matrices
  Z <- list()
  ## for each column of the terms specifying the covariance structure
  counter <- 1
  partitionsR <- list()
  for(j in 1:ncol(Z0)){
    ## for all the random terms decided to be part of the main effect
    # baseZnames <- colnames(init[[1]]$Z)
    for(k in (length(init)-meN+1):(length(init))){
      Z1prov <- init[[k]]$Z # main effect matrix
      cn <- colnames(Z1prov)
      if(is.residual){
        vv <- which(Z0[,j] != 0) # which rows belong to the ith environment
        partitionsR[[j]] <- matrix(c(vv[1],vv[length(vv)]),nrow=1)
        Z[[counter]] <- Z1prov[vv,vv]
        # Z[[counter]] <- Z1prov %*% Diagonal(x=Z0[,j])
      }else{
        provZ0iCol <- Matrix(Z0[,j]) %*% Matrix(1,1,ncol(Z1prov))
        Z1provZ0iCol <- Z1prov * provZ0iCol
        if(!is(class(Z1provZ0iCol), "dgCMatrix")){
          Z[[counter]] <- as(Z1prov * provZ0iCol, Class = "dgCMatrix")
        }else{
          Z[[counter]] <- Z1provZ0iCol
        }
        if(is.null(cn)){
          stop("Gu matrix needs to have row and column names matching the levels of the random effect. Please correct \n", call. = TRUE )
        }
        # print(cn)
        # print(colnames(Gu))
        checkg <- setdiff(cn,colnames(Gu)) # make sure missing levels is not a thing
        # print(checkg)
        if(length(checkg)>0){
          stop(paste("levels of",ref_name,"missing in Gu"),call. = FALSE)
        }
        checkg2 <- setdiff(colnames(Gu),cn) # check if additional G levels exist
        if(length(checkg2)>0){
          cat(paste0("Adding additional levels of Gu in the model matrix of '",ref_name,"' \n"))
          added <- Matrix(0, nrow = nrow(Z1provZ0iCol), ncol = length(checkg2)); colnames(added) <- checkg2
          Z[[counter]] <- cbind(Z1provZ0iCol,added)
        }
      }
      counter <- counter+1
      if(meN <= 1){ # if there's only one main effect
        if(pasteNames & !is.residual){ # if there's a covariance structure paste names, if only isc() don't
          # colnames(Z[[j]]) <- paste(colnames(Z1prov),colnames(Z0)[j],sep=":")
        }
      }
    }
  }
  Zind <- rep(1,length(Z))

  ######################################
  ## meN adjustment
  ## modify theta and thetaC according to the number of mainEffect matrices
  if(is.null(meTheta)){
    meTheta <- diag(meN)
  }
  if(is.null(meThetaC)){
    meThetaC <- diag(meN); rownames(meThetaC) <- colnames(meThetaC) <- letters[1:nrow(meThetaC)]
  }else{meThetaC[lower.tri(meThetaC)]=0}
  if(ncol(meThetaC) > 1){ # if we did a complex structure between effects
    theta <- kronecker(theta,meTheta, make.dimnames = TRUE)
    thetaC <- kronecker(thetaC,meThetaC, make.dimnames = TRUE)
  }

  #########################################
  ## thetaF
  nn <- length(which(thetaC > 0))#unlist(lapply(thetaC, function(x){length(which(x > 0))}))
  # nn2 <- sum(nn[1:max(Zind)])
  thetaF <- diag(nn)
  # print(nrow(thetaF))
  sp0 <- ifelse(sp,1,0)
  sp0 <- rep(sp0,nrow(thetaF))
  if(sp){thetaF <- thetaF*0}
  # we make sure that the A matrix is properly ordered
  cn <- colnames(Z[[length(Z)]])
  Gu <- Gu[cn,cn]
  output <- list(Z=Z, Gu=Gu, theta=theta, thetaC=thetaC, thetaF=thetaF,partitionsR=partitionsR, sp=sp0)
  return(output)
}


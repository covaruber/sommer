vs <- function(..., Gu=NULL, Gti=NULL, Gtc=NULL, reorderGu=TRUE, buildGu=TRUE){
  
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
  init2 <- list() # store the matrices for each random effect provided in ...
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
  if(nre==1){ # if the vs() has only one random effect
    allzs <- matrix(1,nrow=nrow(Z),ncol=1); colnames(allzs) <- "u"
    # if(!is.null(Gtc)){
    #   vcs <- Gtc; colnames(vcs) <- rownames(vcs) <- "u"
    # }else{
    vcs <- matrix(1,1,1); colnames(vcs) <- rownames(vcs) <- "u"
    # }
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
  # print(vcs)
  for(i in 1:ncol(vcs)){ ## for each row
    for(j in 1:i){ ## for each column
      # print(paste(i,j))
      if(vcs[i,j] > 0){ ## to be estimated
        
        if(i==j){## variance component (diagonal term in vcs) and either positive or fixed
          namz <- strsplit(rownames(vcs)[i],":")[[1]]
          zz <- as.matrix(apply(as.matrix(allzs[,namz]),1,prod) * Z)
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

gvs <- function(..., Gu=NULL, Guc=NULL, Gti=NULL, Gtc=NULL, form=NULL){ # general variance structures
  
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
  # ref_name <- namess2[length(namess2)]
  # certain random effects coming from spl2D(), leg(), and others may need some help to find the terms
  specialVariables <- unlist(lapply(init,function(x){(attributes(x)$variables)}))
  # print(namess2)
  if("units" %in% namess2){
    is.residual =TRUE
  }else{is.residual=FALSE}
  ### get the data
  init2 <- list() # store the matrices for each random effect provided in ...
  for(i in 1:length(init)){# i=1
    if(is.list(init[[i]])){ ## if it comes in a list form already from a ds, us, cs function
      
      init2[[i]] <- init[[i]]
      
    }else{ # is a single vector with numbers or characters, ...
      
      if(is.matrix(init[[i]])){ # a matrix is provided already so no need to create it
        
        # relationship matrix
        if(!is.null(Gu)){
          if(length(Gu[[i]]) > 0){
            mm <- Gu[[i]]
          }else{
            mm=diag(ncol(init[[i]])); rownames(mm) <- colnames(mm) <- colnames(init[[i]]) 
          }
        }else{
          mm=diag(ncol(init[[i]])); rownames(mm) <- colnames(mm) <- colnames(init[[i]]) 
        }
        init2[[i]] <- list(x=init[[i]],mm) # store Z and K
        
      }else{ # is a vector
        
        dummy <- init[[i]]
        if(!is.character(dummy) & !is.factor(dummy)){ # user provides a numeric matrix
          dummy <- matrix(dummy,ncol=1) # put the vector in a matrix
          colnames(dummy) <- namess2[i] # add a name to the matrix
          
          if(!is.null(Gu)){
            if(length(Gu[[i]]) > 0){
              mm <- Gu[[i]]
            }else{
              mm=diag(1); rownames(mm) <- colnames(mm) <- namess2[i] # K is a 1 x 1 matrix 
            }
          }else{
            mm=diag(1); rownames(mm) <- colnames(mm) <- namess2[i] # K is a 1 x 1 matrix
          }
          
        }else{ # user provided a factor or character vector, normal case
          levs <- na.omit(unique(dummy)) # extract all levels
          if(length(levs) > 1){ # if more than one level build a design matrix
            dummy  <- model.matrix(~dummy-1,na.action = na.pass) # form Z
          }else{ # there was only one level? wrong but
            vv <- which(!is.na(dummy)); 
            dummy <- matrix(0,nrow=length(dummy)) # create Z
            dummy[vv,] <- 1; colnames(dummy) <- levs # add 1's
          }
          colnames(dummy) <- gsub("dummy","",colnames(dummy)) # add column names
          
          if(!is.null(Gu)){
            if(length(Gu[[i]]) > 0){
              mm <- Gu[[i]]
            }else{
              mm=diag(ncol(dummy)); rownames(mm) <- colnames(mm) <- colnames(dummy) # create K
            }
          }else{
            mm=diag(ncol(dummy)); rownames(mm) <- colnames(mm) <- colnames(dummy) # create K
          }
          
        }
        init2[[i]] <- list(dummy,mm) # store Z and K
      }
    }
  }
  # make a dataframe with the vectors and matrices provided by the user
  # print(str(init2))
  #########################################################
  ## build new design and covariance matrices for the structures
  ## specified in the form argument
  #########################################################
  init3 <- list() # list to store the new matrices
  if(is.null(Guc)){
    Guc <- unsm(length(namess2))
    colnames(Guc) <- rownames(Guc) <- namess2 
  }else{
    colnames(Guc) <- rownames(Guc) <- namess2 
  }
  counter <- 0
  typevc <- numeric()
  re_name <- character()
  for(i in 1:nrow(Guc)){
    for(j in i:ncol(Guc)){
      if(Guc[i,j] != 0){ # if vcov component has to be estimated
        
        counter <- counter+1
        typevc[counter] <- Guc[i,j]
        re_name[counter] <- paste(rownames(Guc)[i],colnames(Guc)[j],sep=":")
        if(i==j){ # variance component
          init3[[counter]] <- list(x=init2[[i]][[1]],init2[[i]][[2]]) # store Z and K
        }else{ # covariance componenent
          Zx <- cbind(init2[[i]][[1]],init2[[j]][[1]])
          Ai <- init2[[i]][[2]]
          Aj <- init2[[j]][[2]]
          Zeroi <- matrix(0, nrow=nrow(Aj),ncol=ncol(Ai)); rownames(Zeroi) <- rownames(Aj); colnames(Zeroi) <- colnames(Ai)
          Zeroj <- matrix(0, nrow=nrow(Ai),ncol=ncol(Aj)); rownames(Zeroj) <- rownames(Ai); colnames(Zeroj) <- colnames(Aj)
          # print(dim(Ai))
          # print(dim(Aj))
          # print(dim(Zeroi))
          # print(dim(Zeroj))
          Ax <- rbind(cbind(Zeroj,Ai),cbind(Aj,Zeroi))
          init3[[counter]] <- list(x=Zx,Ax) # store Z and K 
        }
        
      }
    }
  }
  
  #########################################################
  ## add Gtc and Gti
  #########################################################
  
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
  if(!is.null(specialVariables)){
    namess2 <- specialVariables
  }
  
  Zup <- lapply(init3,function(x){x[[1]]})
  Kup <- lapply(init3,function(x){x[[2]]})
  #########################################################
  ## put it inside a form in case the user wants something 
  ## lie GxE
  #########################################################
  # init4 <- list()
  # if(!is.null(form)){
  #   for(i in 1:length(init3)){
  #     init4[[i]] <- with(data, vs(form, Zup[[i]], Gu=Kup[[i]], Gti=Gti, Gtc=Gtc, reorderGu = FALSE))
  #   }
  #   Zup <- lapply(init4, function(x){x$Z})
  #   Kup <- lapply(init4, function(x){x$K})
  # }
  
  #########################################################
  ## done
  #########################################################
  # Zup <- lapply(init4,function(x){x[[1]]})
  # Kup <- lapply(init4,function(x){x[[2]]})
  
  S3 <- list(Z=Zup,K=Kup,Gti=Gti,Gtc=Gtc,typevc=typevc,re_name=re_name,vcs=Guc, terms=namess2)
  return(S3)
}
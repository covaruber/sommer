#### =========== ######
## PREDICT FUNCTION #
#### =========== ######
# include is used for aggregating
# averaged is used to be included in the prediction
# ignored is not used included in the prediction


"predict.mmes" <- function(object, Dtable=NULL, D, ...){
  if(is.character(D)){classify <- D}else{classify="id"} # save a copy before D is overwriten
  # complete the Dtable withnumber of effects in each term
  xEffectN <- lapply(object$partitionsX, as.vector)
  nz <- unlist(lapply(object$uList,function(x){nrow(x)*ncol(x)}))
  # add a value but if there's no intercept consider it
  lenEff <- length(xEffectN)
  toAdd <- xEffectN[[lenEff]];
  if(length(toAdd) > 0){ # there's a value in the last element of xEffectN
    add <- max(toAdd) + 1 # this is our last column of fixed effects
  }else{ # there's not a value in the last element of xEffectN
    if(lenEff > 1){
      toAdd2 <- xEffectN[[lenEff-1]]
      add <- max(toAdd2) + 1
    }
  }
  if(!is.null(nz)){
    zEffectsN <- list()
    for(i in 1:length(nz)){
      end= add + nz[i] - 1
      zEffectsN[[i]] <- add:end
      add = end + 1
    }
    names(zEffectsN) <- names(nz)
    effectsN = c(xEffectN,zEffectsN)
  }else{
    effectsN = xEffectN
  }
  
  # fill the Dt table for rules
  if(is.null(Dtable) & is.character(D) ){ # if user didn't provide the Dtable but D is character
    Dtable <- object$Dtable # we extract it from the model
    termsInDtable <- apply(data.frame(Dtable$term),1,function(xx){all.names(as.formula(paste0("~",xx)))})
    termsInDtable <- lapply(termsInDtable, function(x){intersect(x,colnames(object$data))})
    termsInDtable <- lapply(termsInDtable, function(x){return(unique(c(x, paste(x, collapse=":"))))})
    # term identified
    termsInDtableN <- unlist(lapply(termsInDtable,length))
    pickTerm <- which( unlist(lapply(termsInDtable, function(xxx){ifelse(length(which(xxx == D)) > 0, 1, 0)})) > 0)
    if(length(pickTerm) == 0){
      isInDf <- which(colnames(object$data) %in% D)
      if(length(isInDf) > 0){ # is in data frame but not in model
        stop(paste("Predict:",classify,"not in the model but present in the original dataset. You may need to provide the
                   Dtable argument to know how to predict", classify), call. = FALSE)
      }else{
        stop(paste("Predict:",classify,"not in the model and not present in the original dataset. Please correct D."), call. = FALSE)
      }
    }
    # check if the term to predict is fixed or random
    pickTermIsFixed = ifelse("fixed" %in% Dtable[pickTerm,"type"], TRUE,FALSE)
    ## 1) when we predict a random effect, fixed effects are purely "average"
    if(!pickTermIsFixed){Dtable[which(Dtable$type %in% "fixed"),"average"]=TRUE; Dtable[pickTerm,"include"]=TRUE}
    ## 2) when we predict a fixed effect, random effects are ignored and the fixed effect is purely "include"
    if(pickTermIsFixed){Dtable[pickTerm,"include"]=TRUE}
    ## 3) for predicting a random effect, the interactions are ignored and only main effect is "included", then we follow 1)
    ## 4) for a model with pure interaction trying to predict a main effect of the interaction we "include" and "average" the interaction and follow 1)
    if(length(pickTerm) == 1){ # only one effect identified
      if(termsInDtableN[pickTerm] > 1){ # we are in #4 (is a pure interaction model)
        Dtable[pickTerm,"average"]=TRUE
      }
    }else{# more than 1, there's main effect and interactions, situation #3
      main <- which(termsInDtableN[pickTerm] == min(termsInDtableN[pickTerm])[1])
      Dtable[pickTerm,"include"]=FALSE;  Dtable[pickTerm,"average"]=FALSE # reset
      Dtable[pickTerm[main],"include"]=TRUE
    }
  }
  ## if user has provided D as a classify then we create the D matrix
  if(is.character(D)){
    # create model matrices to form D
    P <- sparse.model.matrix(as.formula(paste0("~",D,"-1")), data=object$data)
    colnames(P) <- gsub(D,"",colnames(P))
    tP <- t(P)
    W <- object$W
    D = tP %*% W
    colnames(D) <- c(rownames(object$b),rownames(object$u))
    rd <- rownames(D)
    cd <- colnames(D)
    for(jRow in 1:nrow(D)){ # for each effect add 1's where missing
      myMatch <- which(cd == rd[jRow])
      if(length(myMatch) > 0){D[jRow,myMatch]=1}
    }
    # apply rules in Dtable
    for(iRow in 1:nrow(Dtable)){
      w <- effectsN[[iRow]]
      # include/exclude rule
      if(Dtable[iRow,"include"]){ # set to 1
        subD <- D[,w,drop=FALSE]
        subD[which(subD > 0, arr.ind = TRUE)] = 1
        D[,w] <- subD
        # average rule
        if(Dtable[iRow,"average"]){ # set to 1
          # average the include set
          for(o in 1:nrow(subD)){
            v <- which(subD[o,] > 0);  subD[o,v] <- subD[o,v]/length(v)
          }
          D[,w] <- subD
        }
      }else{ # set to zero
        if(Dtable[iRow,"average"]){ # set to 1
          subD <- D[,w,drop=FALSE] + 1
          subD <- subD/subD
          subD[which(subD > 0, arr.ind = TRUE)] = subD[which(subD > 0, arr.ind = TRUE)]/ncol(subD)
          D[,w] <- subD
        }else{
          D[,w] <- D[,w] * 0
        }
      }
    }
    interceptColumn <- unique(c(grep("Intercept",rownames(object$b) ))) # ,which(rownames(object$b)=="1")
    if(length(interceptColumn) > 0){D[,interceptColumn] = 1}
  }else{ }# user has provided D as a matrix to do direct multiplication
  ## calculate predictions and standard errors
  bu <- object$bu
  predicted.value <- D %*% bu
  vcov <- D %*% object$Ci %*% t(D)
  std.error <- sqrt(diag(vcov))
  pvals <- data.frame(id=rownames(D),predicted.value=predicted.value[,1], std.error=std.error)
  if(is.character(classify)){colnames(pvals)[1] <- classify}
  return(list(pvals=pvals,D=D,vcov=vcov, Dtable=Dtable))
}

"print.predict.mmes"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("
                 The predictions are obtained by averaging/aggregating across
                 the hypertable calculated from model terms constructed solely
                 from factors in the include sets. You can customize the model
                 terms used with the 'Dtable' argument.\n")
  ))
  cat(blue(paste("\n Head of predictions:\n")
  ))
  head(x$pvals,...)
}

# "predict.mmes" <- function(object, Dtable=NULL, D, ...){
#
#
#   if(is(D,"character")){ ## if user don't provide a D but a Dtable
#     ## get the names of terms used for each effect from the Dtable
#     termNames <- apply(data.frame(Dtable[,"term"]),1,function(x){
#       y<-paste(all.vars(as.formula(paste("~",x))),collapse = ":")
#       z<-intersect(strsplit(y, split=":")[[1]], colnames(object$data))
#       # z <- z[1] # just the last term to avoid that in interactions more terms are used
#       return(z)
#     })
#     ## get from the Dtable which term will be used for the hypertable classification
#     forD <- D#Dtable[which(Dtable$D),"term"]
#     forD <- paste(all.vars(as.formula(paste("~",forD))),collapse = ":")
#     classify <- intersect(strsplit(forD, split=":")[[1]], colnames(object$data))
#     classifyL <- as.list(classify)
#     object$data[,"Intercept"] <- "1"
#     object$data[,"1"] <- "Intercept"
#     levelsClassifyL <- lapply(classifyL, function(x){sort(unique(object$data[,x]))})
#     ## build a dataset for building the D matrix
#     dataPredict <- do.call(expand.grid, levelsClassifyL )
#     colnames(dataPredict) <- classify
#     ## D matrix
#     Dformed <- Matrix(0, nrow(dataPredict), nrow(object$bu))
#     colnames(Dformed) <- c(rownames(object$b),rownames(object$u))
#
#     #######################
#     ## Fill the fixed part of D
#     #######################
#     DtableF <- Dtable[which(Dtable[,"type"] == "fixed"),]
#     for(i in 1:length(object$partitionsX)){ # for each X partition
#       if( DtableF[i,"include"] ){ # if we want to include
#         for(j in 1:ncol(dataPredict)){
#           for(k in 1:nrow(dataPredict)){
#
#             namesPartition <- rownames(object$b)[object$partitionsX[[i]]]
#             namesPartitionList <- strsplit(namesPartition,split = ":")
#             v <- which( unlist( lapply(namesPartitionList,function(h){
#               res <- ifelse(length(which(h == as.character(dataPredict[k,j]))) > 0,TRUE,FALSE)
#               return(res)
#             })  ))
#
#             if(length(v) > 0){
#               if(DtableF[i,"include"] & !DtableF[i,"average"] ){ # if we want purely include
#                 Dformed[k,namesPartition[v]]=1
#               }
#               if(DtableF[i,"include"] & DtableF[i,"average"]){ # if we want to include and average
#                 Dformed[k,namesPartition[v]]=1/length(v)
#               }
#             }
#
#           }
#         }
#       }else{ # if we don't want to include
#         if(DtableF[i,"average"]){ # but we want to average
#           Dformed[,object$partitionsX[[i]][1,]]=1/length(object$partitionsX[[i]][1,])
#         }
#       }
#     }
#
#     vInt <- which(colnames(Dformed) == "Intercept")
#     if(length(vInt) > 0 ){Dformed[,"Intercept"]=1}
#
#     #######################
#     ## Fill the random part of D
#     #######################
#     DtableR <- Dtable[which(Dtable[,"type"] == "random"),]
#     termNamesR <- termNames[which(Dtable[,"type"] == "random")]
#     for(i in 1:nrow(dataPredict)){ # for each row in the D table or the DataPredict dataset
#
#       dp0 <- as.data.frame(dataPredict[i,]); colnames(dp0) <- colnames(dataPredict)
#       for(j in 1:length(object$partitions)){
#
#         if(DtableR[j,"include"]){
#           # print(dp0)
#
#           dp <- as.data.frame(dp0[,intersect(colnames(dp0),termNamesR[[j]])]); colnames(dp) <-intersect(colnames(dp0),termNamesR[[j]])
#           # print(dp)
#           part <- object$partitions[[j]]
#           partv <- part[1:1]:part[nrow(part),2]
#
#           if(ncol(object$uList[[j]]) == 1){
#             nam2 <- rownames(object$uList[[j]])
#             nam3 <- unlist(lapply(as.list(rownames(object$uList[[j]])),function(x){x2<-strsplit(x,split = ":")[[1]];x2<-x2[length(x2)];return(x2)}))
#           }else{
#             nam2 <- vector(mode="character")
#             for(l in 1:ncol(object$uList[[j]])){
#               nam2 <- c(nam2,paste(colnames(object$uList[[j]])[l],rownames(object$uList[[j]]), sep=":"))
#             }
#             nam3 <- vector(mode="character")
#           }
#
#           namp <- apply( dp , 1 , paste , collapse = ":" )
#           v1 <- which(nam2 %in% namp)
#           v2 <- which(nam3 %in% namp)
#           v <- sort(unique(c(v1,v2)), decreasing = FALSE)
#           if(length(v) > 0){
#             if(DtableR[j,"include"] & !DtableR[j,"average"] ){ # if we want purely include
#               Dformed[i,partv[v]]=1  # Dformed[k,namesPartition[v]]=1
#             }
#             if(DtableR[j,"include"] & DtableR[j,"average"]){ # if we want to include and average
#               Dformed[i,partv[v]]=1/(length(v)*nrow(part)) #Dformed[k,namesPartition[v]]=1/length(v)
#             }
#           }
#         }else{
#           if(DtableR[j,"average"]){ # but we want to average
#             Dformed[,partv]=1/length(partv)
#           }
#         }
#
#
#       }
#     }
#
#   }else if(is(D,"dgCMatrix")){
#     dataPredict <- data.frame(id=1:nrow(D))
#     Dformed <- D
#   }
#
#   Ci <- object$Ci
#
#   vcov <- Dformed%*%Ci%*%t(Dformed)
#   se <- as.vector(sqrt(diag(vcov)))
#
#   predicted.value <- as.vector(Dformed %*% object$bu)
#
#   dataPredict$predicted.value <- predicted.value
#   dataPredict$se <- se
#
#   return(list(pvals=dataPredict, vcov=vcov, D=Dformed, Ci=Ci))
#
# }


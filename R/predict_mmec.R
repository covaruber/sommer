#### =========== ######
## PREDICT FUNCTION #
#### =========== ######
# include is used for aggregating
# averaged is used to be included in the prediction
# ignored is not used included in the prediction

"predict.mmec" <- function(object, Dtable=NULL, D, ...){
  
  
  if(is(D,"character")){ ## if user don't provide a D but a Dtable
    ## get the names of terms used for each effect from the Dtable
    termNames <- apply(data.frame(Dtable[,"term"]),1,function(x){
      y<-paste(all.vars(as.formula(paste("~",x))),collapse = ":")
      z<-intersect(strsplit(y, split=":")[[1]], colnames(object$data))
      # z <- z[1] # just the last term to avoid that in interactions more terms are used
      return(z)
    })
    ## get from the Dtable which term will be used for the hypertable classification
    forD <- D#Dtable[which(Dtable$D),"term"]
    forD <- paste(all.vars(as.formula(paste("~",forD))),collapse = ":")
    classify <- intersect(strsplit(forD, split=":")[[1]], colnames(object$data))
    classifyL <- as.list(classify)
    object$data[,"Intercept"] <- "1"
    object$data[,"1"] <- "Intercept"
    levelsClassifyL <- lapply(classifyL, function(x){sort(unique(object$data[,x]))})
    ## build a dataset for building the D matrix
    dataPredict <- do.call(expand.grid, levelsClassifyL )
    colnames(dataPredict) <- classify
    ## D matrix
    Dformed <- Matrix(0, nrow(dataPredict), nrow(object$bu))
    colnames(Dformed) <- c(rownames(object$b),rownames(object$u))
    
    #######################
    ## Fill the fixed part of D
    #######################
    DtableF <- Dtable[which(Dtable[,"type"] == "fixed"),]
    for(i in 1:length(object$partitionsX)){ # for each X partition
      if( DtableF[i,"include"] ){ # if we want to include
        for(j in 1:ncol(dataPredict)){
          for(k in 1:nrow(dataPredict)){
            
            namesPartition <- rownames(object$b)[object$partitionsX[[i]]]
            namesPartitionList <- strsplit(namesPartition,split = ":")
            v <- which( unlist( lapply(namesPartitionList,function(h){
              res <- ifelse(length(which(h == as.character(dataPredict[k,j]))) > 0,TRUE,FALSE)
              return(res)
            })  ))
            
            if(length(v) > 0){
              if(DtableF[i,"include"] & !DtableF[i,"average"] ){ # if we want purely include
                Dformed[k,namesPartition[v]]=1 
              }
              if(DtableF[i,"include"] & DtableF[i,"average"]){ # if we want to include and average
                Dformed[k,namesPartition[v]]=1/length(v)
              }
            }
            
          }
        }
      }else{ # if we don't want to include
        if(DtableF[i,"average"]){ # but we want to average
          Dformed[,object$partitionsX[[i]][1,]]=1/length(object$partitionsX[[i]][1,])
        }
      }
    }
    
    vInt <- which(colnames(Dformed) == "Intercept")
    if(length(vInt) > 0 ){Dformed[,"Intercept"]=1}
    
    #######################
    ## Fill the random part of D
    #######################
    DtableR <- Dtable[which(Dtable[,"type"] == "random"),]
    termNamesR <- termNames[which(Dtable[,"type"] == "random")]
    for(i in 1:nrow(dataPredict)){ # for each row in the D table or the DataPredict dataset
      
      dp0 <- as.data.frame(dataPredict[i,]); colnames(dp0) <- colnames(dataPredict)
      for(j in 1:length(object$partitions)){
        
        if(DtableR[j,"include"]){
          # print(dp0)
          
          dp <- as.data.frame(dp0[,intersect(colnames(dp0),termNamesR[[j]])]); colnames(dp) <-intersect(colnames(dp0),termNamesR[[j]])
          # print(dp)
          part <- object$partitions[[j]]
          partv <- part[1:1]:part[nrow(part),2]
          
          if(ncol(object$uList[[j]]) == 1){
            nam2 <- rownames(object$uList[[j]])
            nam3 <- unlist(lapply(as.list(rownames(object$uList[[j]])),function(x){x2<-strsplit(x,split = ":")[[1]];x2<-x2[length(x2)];return(x2)}))
          }else{
            nam2 <- vector(mode="character")
            for(l in 1:ncol(object$uList[[j]])){
              nam2 <- c(nam2,paste(colnames(object$uList[[j]])[l],rownames(object$uList[[j]]), sep=":"))
            }
            nam3 <- vector(mode="character")
          }
          
          namp <- apply( dp , 1 , paste , collapse = ":" )
          v1 <- which(nam2 %in% namp)
          v2 <- which(nam3 %in% namp)
          v <- sort(unique(c(v1,v2)), decreasing = FALSE)
          if(length(v) > 0){
            if(DtableR[j,"include"] & !DtableR[j,"average"] ){ # if we want purely include
              Dformed[i,partv[v]]=1  # Dformed[k,namesPartition[v]]=1 
            }
            if(DtableR[j,"include"] & DtableR[j,"average"]){ # if we want to include and average
              Dformed[i,partv[v]]=1/(length(v)*nrow(part)) #Dformed[k,namesPartition[v]]=1/length(v)
            }
          }
        }else{
          if(DtableR[j,"average"]){ # but we want to average
            Dformed[,partv]=1/length(partv)
          }
        }
        
        
      }
    }
    
  }else if(is(D,"dgCMatrix")){ 
    dataPredict <- data.frame(id=1:nrow(D))
    Dformed <- D
  }
  
  Ci <- object$Ci
  
  vcov <- Dformed%*%Ci%*%t(Dformed)
  se <- as.vector(sqrt(diag(vcov))) 
  
  predicted.value <- as.vector(Dformed %*% object$bu)
  
  dataPredict$predicted.value <- predicted.value
  dataPredict$se <- se
  
  return(list(pvals=dataPredict, vcov=vcov, D=Dformed, Ci=Ci))
  
}

"print.predict.mmec"<- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(blue(paste("
                 The predictions are obtained by averaging/aggregating across
                 the hypertable calculated from model terms constructed solely
                 from factors in the include sets. You can customize the model
                 terms used with the 'hypertable' argument. Current model terms used:\n")
  ))
  cat(blue(paste("\n Head of predictions:\n")
  ))
  head(x$pvals,...)
}